"""Poisson placement and process-based parallel family simulation."""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from itertools import islice
import math
import multiprocessing as mp
import numpy as np

from .certificate import Certificate, certify
from .model import Model, finite, integer
from .trees import Family, Skeleton, ResourceLimitError, make_family, make_skeleton, rng_for, _readonly


@dataclass(frozen=True)
class Plan:
    model: Model
    certificate: Certificate
    horizon: float
    seed: int
    arrivals: np.ndarray


@dataclass(frozen=True)
class Forest:
    """Natural FRIME lifetimes on [0,horizon], initialised using a finite past.

    T truncates immigration HISTORY, never a family's lifetime. `at` is
    deterministic: queries at different times share the same realised forest.
    """
    plan: Plan
    families: tuple[Family, ...]

    def at(self, t: float = 0.0) -> np.ndarray:
        finite(t, "t")
        if not 0 <= t <= self.plan.horizon:
            raise ValueError("t must lie in [0, horizon]")
        parts = [family.at(t-float(s)) for s, family in zip(self.plan.arrivals, self.families) if s <= t]
        return np.concatenate(parts) if parts else np.empty(0)

    def histogram(self, bins, t: float = 0.0) -> np.ndarray:
        """Bin raw fragment COUNTS, not a normalised density."""
        finite(t, "t")
        if not 0 <= t <= self.plan.horizon:
            raise ValueError("t must lie in [0, horizon]")
        edges = np.asarray(bins, dtype=float)
        if edges.ndim != 1 or len(edges) < 2 or not np.all(np.isfinite(edges)) or not np.all(np.diff(edges) > 0):
            raise ValueError("bins must be finite and strictly increasing")
        result = np.zeros(len(edges)-1, dtype=np.int64)
        for s, family in zip(self.plan.arrivals, self.families):
            if s <= t:
                result += np.histogram(family.at(t-float(s)), edges)[0]
        return result

    @property
    def stored_nodes(self) -> int:
        return sum(len(f.sizes) for f in self.families)


def make_plan(model: Model, *, delta: float = 1e-6, horizon: float = 0.0, seed: int = 0,
              max_families: int = 100_000) -> Plan:
    """Draw a Poisson process on [-T,horizon]. Unsorted times are sufficient."""
    finite(horizon, "horizon")
    if horizon < 0:
        raise ValueError("horizon must be nonnegative")
    seed = integer(seed, "seed")
    max_families = integer(max_families, "max_families", 1)
    certificate = certify(model, delta)
    width = certificate.T + horizon
    mean = model.c_i * width
    if not math.isfinite(width) or not math.isfinite(mean) or mean > max_families:
        raise ResourceLimitError(f"expected immigrant count {mean:g} exceeds max_families={max_families}")
    rng = rng_for(seed, 0, 3)
    n = int(rng.poisson(mean))
    if n > max_families:
        raise ResourceLimitError(f"realised immigrant count {n} exceeds max_families={max_families}")
    arrivals = rng.uniform(-certificate.T, horizon, n) if n else np.empty(0)
    return Plan(model, certificate, horizon, seed, _readonly(arrivals))


def _batch(task):
    model, seed, jobs, mode, max_nodes = task
    if mode == "prune":
        return [tree.prune(model) for _, tree in jobs]
    if mode in {"skeleton", "marked"}:
        trees = [make_skeleton(model, i, seed, max_nodes=max_nodes) for i, _ in jobs]
        return [tree.prune(model) for tree in trees] if mode == "marked" else trees
    return [make_family(model, i, seed, age, max_nodes=max_nodes) for i, age in jobs]


def _map_batches(model, seed, jobs, mode, workers, batch_size, max_nodes, max_total_nodes):
    workers = integer(workers, "workers", 1)
    batch_size = integer(batch_size, "batch_size", 1)
    max_nodes = integer(max_nodes, "max_nodes", 1)
    max_total_nodes = integer(max_total_nodes, "max_total_nodes", 1)
    iterator = iter(jobs)
    def tasks():
        while chunk := list(islice(iterator, batch_size)):
            yield (model, seed, chunk, mode, max_nodes)
    output, stored = [], 0
    def append_batch(batch):
        nonlocal stored
        stored += sum(len(f.sizes) for f in batch)
        if stored > max_total_nodes:
            raise ResourceLimitError(f"forest exceeds max_total_nodes={max_total_nodes}")
        output.extend(batch)
    if workers == 1:
        for task in tasks():
            append_batch(_batch(task))
    else:
        # Explicit spawn works on macOS and avoids reliance on a forked RNG state.
        # Keep only workers batches in flight, bounding queued/returned data.
        with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as pool:
            source = iter(tasks())
            pending = [pool.submit(_batch, task) for task in islice(source, workers)]
            while pending:
                append_batch(pending.pop(0).result())
                task = next(source, None)
                if task is not None:
                    pending.append(pool.submit(_batch, task))
    return tuple(output)


def presample(model: Model, n: int, *, seed: int = 0, workers: int = 1, batch_size: int = 32,
              max_nodes: int = 1_000_000, max_total_nodes: int = 2_000_000) -> tuple[Skeleton, ...]:
    """Generate n independent pure trees. Each tree is used once per population.

    A bank can be reused for paired parameter comparisons, but those resulting
    populations are correlated and are NOT independent Monte Carlo replicates.
    """
    n = integer(n, "n")
    integer(seed, "seed")
    if n > 100_000:
        raise ResourceLimitError("presample supports at most 100000 trees per bank")
    return _map_batches(model, seed, ((i, 0.0) for i in range(n)), "skeleton", workers,
                        batch_size, max_nodes, max_total_nodes)


def assemble(plan: Plan, bank: tuple[Skeleton, ...], *, workers: int = 1, batch_size: int = 32,
             max_total_nodes: int = 2_000_000) -> Forest:
    """Attach one fresh independent bank entry to each Poisson point.

    Caller must use an iid bank independent of plan.arrivals; do not select
    trees by survival, size, or results. make_plan/presample use disjoint streams.
    """
    if len(bank) < len(plan.arrivals):
        raise ValueError("bank exhausted: generate more independent trees, do not recycle")
    chosen = bank[:len(plan.arrivals)]
    if len({s.family_id for s in chosen}) != len(chosen):
        raise ValueError("duplicate family IDs would reuse randomness within a population")
    families = _map_batches(plan.model, 0, enumerate(chosen), "prune", workers,
                            batch_size, 1, max_total_nodes)
    return Forest(plan, families)


def sample(model: Model, *, delta: float = 1e-6, horizon: float = 0.0, seed: int = 0,
           workers: int = 1, batch_size: int = 32, mode: str = "lazy",
           max_families: int = 100_000, max_nodes: int = 1_000_000,
           max_total_nodes: int = 2_000_000) -> Forest:
    """Sample an approximately stationary forest with a reported history bound."""
    if mode not in {"lazy", "skeleton"}:
        raise ValueError("mode must be lazy or skeleton")
    plan = make_plan(model, delta=delta, horizon=horizon, seed=seed, max_families=max_families)
    jobs = ((i, horizon-float(s)) for i, s in enumerate(plan.arrivals))
    execution = "marked" if mode == "skeleton" else "lazy"
    families = _map_batches(model, seed, jobs, execution, workers, batch_size, max_nodes, max_total_nodes)
    return Forest(plan, families)
