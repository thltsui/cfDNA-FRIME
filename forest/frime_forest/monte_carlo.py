"""Independent stationary-population Monte Carlo with persistent family workers."""
from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from datetime import datetime, timezone
from itertools import islice
import math
import multiprocessing as mp
from time import perf_counter

import numpy as np

from .certificate import Certificate
from .model import Model, integer
from .simulation import _batch, make_plan
from .trees import ResourceLimitError, _readonly


def replicate_seed(seed: int, replicate_id: int) -> int:
    """Injectively encode (root seed, global replicate ID) in an MC namespace.

    Cantor pairing avoids overlaps between ensembles with adjacent root seeds.
    The existing SeedSequence/PCG64 machinery mixes this integer into streams.
    Repeating the same pair deliberately reproduces the same population.
    """
    seed = integer(seed, "seed")
    replicate_id = integer(replicate_id, "replicate_id")
    total = seed + replicate_id
    pair = total * (total + 1) // 2 + replicate_id
    return (pair << 32) | 0x4D435F31  # ASCII MC_1: versioned stream namespace


@dataclass(frozen=True)
class MonteCarloResult:
    model: Model
    certificate: Certificate
    seed: int
    replicate_ids: np.ndarray
    seeds: tuple[int, ...]
    bins: np.ndarray
    histograms: np.ndarray
    counts: np.ndarray
    masses: np.ndarray
    immigrants: np.ndarray
    generated_nodes: np.ndarray
    seconds: np.ndarray
    elapsed_seconds: float
    started_utc: str
    workers: int
    batch_size: int
    mode: str

    @property
    def joint_history_bound(self) -> float:
        """Ideal-law mismatch bound for the whole independent replicate vector."""
        return -math.expm1(-len(self.seeds) * self.certificate.lambda_bound)

    def statistics(self) -> dict:
        """Replicate means and Monte Carlo SEs; no SE is available for n=1.

        These estimate finite-history expectations. A history TV bound is not
        an absolute bias bound for unbounded observables such as abundance.
        """
        result = {}
        for name, values in (("histogram", self.histograms), ("count", self.counts),
                             ("mass", self.masses)):
            with np.errstate(over="raise", invalid="raise"):
                mean = np.asarray(values.mean(axis=0)).tolist()
                se = (np.asarray(values.std(axis=0, ddof=1)/math.sqrt(len(values))).tolist()
                      if len(values) > 1 else None)
            result[name] = dict(mean=mean, mcse=se)
        return result


def _summary_batch(task):
    m, seed, jobs, execution, max_nodes, bins, max_total_nodes = task
    families = _batch((m, seed, jobs, execution, max_nodes))
    nodes = sum(len(f.sizes) for f in families)
    if nodes > max_total_nodes:
        raise ResourceLimitError(f"batch exceeds max_total_nodes={max_total_nodes}")
    values = [f.at(age) for f, (_, age) in zip(families, jobs)]
    masses = [float(x.sum()) for x in values]
    if not all(math.isfinite(mass) for mass in masses):
        raise FloatingPointError("retained mass is not representable")
    joined = np.concatenate(values) if values else np.empty(0)
    # Return summaries, not genealogies. Preserve per-family masses so changing
    # batch boundaries does not change floating-point reduction order.
    return np.histogram(joined, bins)[0], len(joined), masses, nodes


class MonteCarlo:
    """Context-managed persistent workers across immigrant-family batches.

    run() uses fresh streams indexed by (seed, start_replicate + row).
    Keep delta fixed when extending an ensemble in chunks. Reusing IDs is
    deterministic replay, not a fresh replicate. Calls on one object are serial.
    """

    def __init__(self, *, workers: int = 1, batch_size: int = 128,
                 max_families: int = 100_000, max_nodes: int = 1_000_000,
                 max_total_nodes: int = 2_000_000):
        self.workers = integer(workers, "workers", 1)
        self.batch_size = integer(batch_size, "batch_size", 1)
        self.max_families = integer(max_families, "max_families", 1)
        self.max_nodes = integer(max_nodes, "max_nodes", 1)
        self.max_total_nodes = integer(max_total_nodes, "max_total_nodes", 1)
        self._pool = None
        self._active = False

    def __enter__(self):
        if self._active:
            raise RuntimeError("MonteCarlo context is already active")
        if self.workers > 1:
            self._pool = ProcessPoolExecutor(max_workers=self.workers,
                                            mp_context=mp.get_context("spawn"))
        self._active = True
        return self

    def close(self):
        if self._pool is not None:
            self._pool.shutdown(wait=True, cancel_futures=True)
            self._pool = None
        self._active = False

    def __exit__(self, *exc):
        self.close()

    def _summarize(self, plan, bins, mode):
        jobs = iter((i, -float(arrival)) for i, arrival in enumerate(plan.arrivals))
        execution = "marked" if mode == "skeleton" else "lazy"
        def tasks():
            while chunk := list(islice(jobs, self.batch_size)):
                yield (plan.model, plan.seed, chunk, execution, self.max_nodes,
                       bins, self.max_total_nodes)
        histogram = np.zeros(len(bins)-1, dtype=np.int64)
        count, nodes, masses = 0, 0, []
        def consume(result):
            nonlocal count, nodes
            h, n, per_family_mass, generated = result
            nodes += generated
            if nodes > self.max_total_nodes:
                raise ResourceLimitError(f"population exceeds max_total_nodes={self.max_total_nodes}")
            histogram[:] += h
            count += n
            masses.extend(per_family_mass)
        if self._pool is None:
            for task in tasks():
                consume(_summary_batch(task))
        else:
            source = iter(tasks())
            pending = [self._pool.submit(_summary_batch, task)
                       for task in islice(source, self.workers)]
            try:
                while pending:
                    consume(pending.pop(0).result())
                    task = next(source, None)
                    if task is not None:
                        pending.append(self._pool.submit(_summary_batch, task))
            finally:
                for future in pending:
                    future.cancel()
        try:
            mass = math.fsum(masses)
        except OverflowError as exc:
            raise FloatingPointError("population mass is not representable") from exc
        return histogram, count, mass, nodes

    def run(self, model: Model, n: int, *, bins, seed: int = 0,
            start_replicate: int = 0, delta: float = 1e-6,
            mode: str = "lazy") -> MonteCarloResult:
        """Generate iid finite-history population summaries at observation time 0.

        delta is PER REPLICATE, not a tolerance for the complete ensemble.
        Caps apply per population and errors abort the whole call without
        returning successful prefixes. No failed seed is skipped or retried.
        Bins may cover a subwindow; counts/masses always cover all retained sizes.
        """
        if not self._active:
            raise RuntimeError("use MonteCarlo inside a with block")
        n = integer(n, "n", 1)
        seed = integer(seed, "seed")
        start_replicate = integer(start_replicate, "start_replicate")
        if start_replicate + n > 2**63:
            raise ValueError("replicate IDs must fit signed int64")
        if mode not in {"lazy", "skeleton"}:
            raise ValueError("mode must be lazy or skeleton")
        edges = np.array(bins, dtype=float, copy=True)
        if (edges.ndim != 1 or len(edges) < 2 or not np.all(np.isfinite(edges))
                or not np.all(np.diff(edges) > 0)):
            raise ValueError("bins must be finite and strictly increasing")
        ids = np.arange(start_replicate, start_replicate+n, dtype=np.int64)
        seeds = tuple(replicate_seed(seed, int(r)) for r in ids)
        histograms = np.empty((n, len(edges)-1), dtype=np.int64)
        counts, immigrants, nodes = (np.empty(n, dtype=np.int64) for _ in range(3))
        masses, seconds = np.empty(n), np.empty(n)
        timestamp = datetime.now(timezone.utc).isoformat()
        start = perf_counter()
        for row, derived_seed in enumerate(seeds):
            tick = perf_counter()
            try:
                plan = make_plan(model, delta=delta, seed=derived_seed,
                                 max_families=self.max_families)
                h, count, mass, generated = self._summarize(plan, edges, mode)
            except (ResourceLimitError, FloatingPointError) as exc:
                self.close()
                raise type(exc)(f"replicate {int(ids[row])}: {exc}") from exc
            except BaseException:
                self.close()
                raise
            histograms[row], counts[row], masses[row] = h, count, mass
            immigrants[row], nodes[row] = len(plan.arrivals), generated
            seconds[row] = perf_counter()-tick
        elapsed = perf_counter()-start
        return MonteCarloResult(model, plan.certificate, seed, _readonly(ids, np.int64), seeds,
                                _readonly(edges), _readonly(histograms, np.int64),
                                _readonly(counts, np.int64), _readonly(masses),
                                _readonly(immigrants, np.int64), _readonly(nodes, np.int64),
                                _readonly(seconds), elapsed, timestamp, self.workers,
                                self.batch_size, mode)


def monte_carlo(model: Model, n: int, *, bins, seed: int = 0,
                start_replicate: int = 0, delta: float = 1e-6, mode: str = "lazy",
                workers: int = 1, batch_size: int = 128,
                max_families: int = 100_000, max_nodes: int = 1_000_000,
                max_total_nodes: int = 2_000_000) -> MonteCarloResult:
    """One ensemble with one pool. Use MonteCarlo's context for repeated calls."""
    with MonteCarlo(workers=workers, batch_size=batch_size, max_families=max_families,
                    max_nodes=max_nodes, max_total_nodes=max_total_nodes) as runner:
        return runner.run(model, n, bins=bins, seed=seed, start_replicate=start_replicate,
                          delta=delta, mode=mode)
