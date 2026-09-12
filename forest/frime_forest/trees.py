"""Independent family marks, pure skeletons, and inherited exit pruning."""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np

from .model import Model, finite, integer


class ResourceLimitError(RuntimeError):
    """A resource cap aborts the run; no silently truncated sample is returned."""


def rng_for(seed: int, family_id: int, component: int = 0) -> np.random.Generator:
    seed = integer(seed, "seed")
    family_id = integer(family_id, "family_id")
    component = integer(component, "component")
    return np.random.Generator(np.random.PCG64(np.random.SeedSequence([component, family_id, seed])))


def _readonly(values, dtype=float):
    result = np.asarray(values, dtype=dtype)
    result.setflags(write=False)
    return result


def _advance(birth: float, wait: float) -> float:
    end = birth + wait
    if not math.isfinite(end) or end <= birth:
        raise FloatingPointError("event time under/overflows; rescale time or length")
    return end


def _children(x: float, r: float):
    left = x*r
    right = x-left
    if not (0 < left < x and 0 < right < x):
        raise FloatingPointError("Beta split rounded to an endpoint; use less extreme parameters")
    return left, right


@dataclass(frozen=True)
class Family:
    family_id: int
    root_size: float
    sizes: np.ndarray
    births: np.ndarray
    ends: np.ndarray
    parents: np.ndarray
    horizon: float

    def at(self, age: float) -> np.ndarray:
        finite(age, "age")
        if age < 0 or age > self.horizon:
            raise ValueError("age must lie in the recorded family horizon")
        return self.sizes[(self.births <= age) & (age < self.ends)].copy()


@dataclass(frozen=True)
class Skeleton:
    """Complete pure-fragmentation tree above ell, with independent exit marks.

    Array order is parent-before-child. Potential descendants remain present
    even when an exit mechanism will later prune an ancestor.
    """
    family_id: int
    root_size: float
    fragmentation_key: tuple
    sizes: np.ndarray
    births: np.ndarray
    splits: np.ndarray
    parents: np.ndarray
    exit_marks: np.ndarray

    def prune(self, model: Model) -> Family:
        if model.fragmentation_key() != self.fragmentation_key:
            raise ValueError("skeleton fragmentation parameters differ; only exit/c_i may change")
        n = len(self.sizes)
        reachable = np.zeros(n, dtype=bool)
        splits_alive = np.zeros(n, dtype=bool)
        ends = np.empty(n)
        for i in range(n):
            parent = self.parents[i]
            reachable[i] = parent < 0 or splits_alive[parent]
            if not reachable[i]:
                ends[i] = self.births[i]
                continue
            rate = model.exit.rate(float(self.sizes[i]), model.L)
            kill = (_advance(float(self.births[i]), float(self.exit_marks[i])/rate)
                    if rate else math.inf)
            ends[i] = min(kill, self.splits[i])
            splits_alive[i] = self.splits[i] < kill
        # Preserve indices so parent IDs remain interpretable, even for pruned nodes.
        return Family(self.family_id, self.root_size, self.sizes, self.births,
                      _readonly(ends), self.parents, math.inf)


def make_skeleton(model: Model, family_id: int, seed: int, *, root_size: float | None = None,
                  max_nodes: int = 1_000_000) -> Skeleton:
    max_nodes = integer(max_nodes, "max_nodes", 1)
    rng = rng_for(seed, family_id, 1)
    x0 = float(rng_for(seed, family_id, 2).uniform(0, model.L)) if root_size is None else finite(root_size, "root_size")
    if not 0 <= x0 <= model.L:
        raise ValueError("root_size must lie in [0,L]")
    sizes, births, splits, parents, marks = [], [], [], [], []
    stack = [(x0, 0.0, -1)] if x0 > model.ell else []
    while stack:
        if len(sizes) >= max_nodes:
            raise ResourceLimitError(f"family {family_id} exceeds max_nodes={max_nodes}")
        x, birth, parent = stack.pop()
        r = float(rng.beta(model.a, model.b))
        split = _advance(birth, float(rng.exponential())/model.rho(x))
        mark = float(rng.exponential())
        if mark <= 0:
            raise FloatingPointError("zero exit mark is not representable in the ideal model")
        index = len(sizes)
        sizes.append(x); births.append(birth); splits.append(split); parents.append(parent); marks.append(mark)
        left, right = _children(x, r)
        if right > model.ell:
            stack.append((right, split, index))
        if left > model.ell:
            stack.append((left, split, index))
    return Skeleton(family_id, x0, model.fragmentation_key(), _readonly(sizes),
                    _readonly(births), _readonly(splits), _readonly(parents, np.int64), _readonly(marks))


def make_family(model: Model, family_id: int, seed: int, horizon: float, *,
                max_nodes: int = 1_000_000) -> Family:
    """Lazy genealogy: stop at exit, size cutoff, or the observation horizon.

    Uses the same family-level stream as skeleton generation. Traversal changes
    after pruning, so the two modes need not give pathwise identical samples.
    """
    finite(horizon, "horizon")
    if horizon < 0:
        raise ValueError("horizon must be nonnegative")
    max_nodes = integer(max_nodes, "max_nodes", 1)
    rng = rng_for(seed, family_id, 1)
    x0 = float(rng_for(seed, family_id, 2).uniform(0, model.L))
    sizes, births, ends, parents = [], [], [], []
    stack = [(x0, 0.0, -1)] if x0 > model.ell else []
    while stack:
        if len(sizes) >= max_nodes:
            raise ResourceLimitError(f"family {family_id} exceeds max_nodes={max_nodes}")
        x, birth, parent = stack.pop()
        r = float(rng.beta(model.a, model.b))
        split = _advance(birth, float(rng.exponential())/model.rho(x))
        mark = float(rng.exponential())
        if mark <= 0:
            raise FloatingPointError("zero exit mark")
        rate = model.exit.rate(x, model.L)
        kill = _advance(birth, mark/rate) if rate else math.inf
        end = min(split, kill)
        index = len(sizes)
        sizes.append(x); births.append(birth); ends.append(end); parents.append(parent)
        if split >= kill or split > horizon:
            continue
        left, right = _children(x, r)
        if right > model.ell:
            stack.append((right, split, index))
        if left > model.ell:
            stack.append((left, split, index))
    return Family(family_id, x0, _readonly(sizes), _readonly(births), _readonly(ends),
                  _readonly(parents, np.int64), horizon)
