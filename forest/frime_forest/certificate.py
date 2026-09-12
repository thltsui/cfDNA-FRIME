"""Non-asymptotic past-history error bound derived in docs/MATHEMATICS.md.

The certificate concerns ideal-law history truncation. It is not a bound on
floating-point or pseudorandom-number error and is not an empirical CI.
"""
from __future__ import annotations

from dataclasses import dataclass
import math
import numpy as np

from .model import Model, finite, integer


@dataclass(frozen=True)
class Certificate:
    model: Model
    delta: float
    T: float
    p: int
    kappa: float
    log_A: float
    lambda_bound: float
    mismatch_bound: float
    intervals: int


def contraction(a: float, b: float, p: int) -> float:
    """1 - E[R**p+(1-R)**p] for integer p>=2, evaluated with log products."""
    integer(p, "p", 2)
    if p == 2:
        s = a + b
        result = 2 * (a / s) * (b / s) * (s / (s + 1))
    else:
        left = math.exp(math.fsum(math.log(a+j)-math.log(a+b+j) for j in range(p)))
        right = math.exp(math.fsum(math.log(b+j)-math.log(a+b+j) for j in range(p)))
        result = 1.0 - math.fsum((left, right))
    if not 0 < result < 1:
        raise FloatingPointError("split contraction is not numerically resolvable")
    return result


def certify(model: Model, delta: float = 1e-6, *, powers=(2, 3, 4, 6, 8),
            intervals: int = 256) -> Certificate:
    """Choose the shortest certified lookback among the specified moment bounds.

    On each size interval, inf(f+g) >= inf(f)+inf(g). Power fragmentation is
    monotone, and supported exit functions are nonincreasing. Endpoint minima
    therefore give a LOWER bound, not a sampled approximation to an infimum.
    """
    finite(delta, "delta")
    if not 0 < delta < 1:
        raise ValueError("delta must lie strictly between zero and one")
    intervals = integer(intervals, "intervals", 1)
    powers = tuple(integer(p, "p", 2) for p in powers)
    if not powers or max(powers) > 100:
        raise ValueError("supply at least one integer moment between 2 and 100")
    edges = np.geomspace(model.ell, model.L, intervals + 1)
    # Ensure exact coverage despite roundoff in generated endpoints.
    edges[0], edges[-1] = model.ell, model.L
    lo, hi = edges[:-1], edges[1:]
    min_rho = np.array([min(model.rho(float(x)), model.rho(float(y))) for x,y in zip(lo,hi)])
    min_exit = np.array([model.exit.rate(float(x), model.L) for x in hi])
    log_ratio = math.log(model.L) - math.log(model.ell)
    candidates = []
    for p in powers:
        kappa = float(np.min(contraction(model.a, model.b, p) * min_rho + min_exit))
        # Small roundoff margin; this is not formal interval arithmetic.
        kappa *= 1 - 1e-12
        if not math.isfinite(kappa) or kappa <= 0:
            raise FloatingPointError("positive decay bound is not representable")
        log_A = p*log_ratio - math.log(p+1) + math.log(-math.expm1(-(p+1)*log_ratio))
        if model.c_i == 0:
            T, lam = 0.0, 0.0
        else:
            log_B = math.log(model.c_i) + log_A - math.log(kappa)
            # Lambda <= delta is slightly stronger than 1-exp(-Lambda) <= delta.
            T = max(0.0, (log_B - math.log(delta))/kappa)
            if not math.isfinite(T):
                raise FloatingPointError("lookback overflows; rescale the model")
            T = float(np.nextafter(T, math.inf))
            lam = math.exp(log_B - kappa*T)
        candidates.append(Certificate(model, delta, T, p, kappa, log_A, lam,
                                      -math.expm1(-lam), intervals))
    return min(candidates, key=lambda c: c.T)
