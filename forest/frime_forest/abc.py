"""Observation and rejection helpers for the exploratory length-profile ABC pilot.

These functions do not replace population simulation with a mean profile.
"""
from __future__ import annotations
import numpy as np
from .model import integer


def counts_array(values):
    x = np.asarray(values)
    if (x.ndim != 1 or len(x) < 1 or not np.issubdtype(x.dtype, np.number)
            or not np.all(np.isfinite(x)) or np.any(x < 0)
            or np.any(x != np.floor(x)) or np.any(x >= 2**53)):
        raise ValueError("counts must be a nonempty vector of finite nonnegative exact integers")
    if float(x.sum(dtype=float)) >= 2**53:
        raise ValueError("total count is too large for exact float64 normalization")
    return x.astype(np.int64)


def observe_counts(population_histogram, n, rng):
    """Equal-detectability multinomial reads; empty population returns no reads.

    The all-zero outcome is a failure category, never a request to resimulate.
    This is a working observation model conditional on a planned read depth.
    """
    h = counts_array(population_histogram)
    n = integer(n, "n", 1)
    if n >= 2**53:
        raise ValueError("read depth is too large")
    return rng.multinomial(n, h / h.sum()) if h.sum() else np.zeros_like(h)


def cdf_distance(simulated, observed):
    """Maximum CDF gap on the supplied ordered bins; no-read outcome is infinity."""
    x, y = counts_array(simulated), counts_array(observed)
    if x.shape != y.shape or y.sum() == 0:
        raise ValueError("matching bins and positive observed read depth are required")
    if x.sum() == 0:
        return float("inf")
    return float(np.max(np.abs(np.cumsum(x) / x.sum() - np.cumsum(y) / y.sum())))


def rejection_indices(distances, epsilon):
    """Keep finite distances <= a frozen, finite tolerance. Never keep no-read outcomes."""
    d = np.asarray(distances, dtype=float)
    if d.ndim != 1 or np.any(np.isnan(d)) or np.any(d < 0):
        raise ValueError("distances must be nonnegative, one-dimensional and not NaN")
    if not np.isfinite(epsilon) or not 0 <= epsilon <= 1:
        raise ValueError("epsilon must be finite and in [0,1]")
    return np.flatnonzero(np.isfinite(d) & (d <= epsilon))
