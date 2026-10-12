"""Thesis Chapter 7 conventions; no dependence on the legacy simulator."""
from __future__ import annotations

from dataclasses import dataclass, field
import math
import numbers


def integer(value: int, name: str, minimum: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, numbers.Integral) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(value)


def finite(value: float, name: str, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, numbers.Real):
        raise ValueError(f"{name} must be a real number")
    value = float(value)
    if not math.isfinite(value) or (positive and value <= 0):
        raise ValueError(f"{name} must be finite" + (" and positive" if positive else ""))
    return value


@dataclass(frozen=True)
class Exit:
    """CON, PNB or PFB. c_e=1 reproduces the thesis power-exit normalisation."""

    kind: str = "constant"
    c_e: float = 1.0
    power: float = -1.0
    boundary: float = 0.4

    def __post_init__(self) -> None:
        if self.kind not in {"constant", "pnb", "pfb"}:
            raise ValueError("exit kind must be constant, pnb or pfb")
        if finite(self.c_e, "c_e") < 0:
            raise ValueError("c_e must be nonnegative")
        if finite(self.power, "exit power") >= 0 and self.kind != "constant":
            raise ValueError("power exits require a negative exponent")
        finite(self.boundary, "boundary", positive=True)

    def rate(self, x: float, L: float) -> float:
        if self.kind == "constant" or self.c_e == 0:
            return float(self.c_e)
        if self.kind == "pfb" and x >= self.boundary:
            return 0.0
        rate = self.c_e * math.exp(self.power * (math.log(x) - math.log(L)))
        if self.kind == "pfb":
            # Stable evaluation of x**power - boundary**power near the boundary.
            rate *= -math.expm1(self.power * (math.log(self.boundary) - math.log(x)))
        # Exact zero exits returned above; a computed zero here is underflow.
        if not math.isfinite(rate) or rate <= 0:
            raise FloatingPointError("exit rate is not representable; rescale the model")
        return rate


@dataclass(frozen=True)
class Immigration:
    """Normalized immigrant root law, separate from split-ratio parameters."""

    kind: str = "uniform"
    a: float = 1.0
    b: float = 1.0

    def __post_init__(self) -> None:
        if self.kind not in {"uniform", "beta"}:
            raise ValueError("immigration kind must be uniform or beta")
        object.__setattr__(self, "a", finite(self.a, "immigration a", positive=True))
        object.__setattr__(self, "b", finite(self.b, "immigration b", positive=True))
        if not math.isfinite(self.a + self.b):
            raise ValueError("immigration a+b must be representable")
        if self.kind == "uniform" and (self.a != 1 or self.b != 1):
            raise ValueError("uniform immigration requires unit shapes")

    @property
    def is_uniform(self) -> bool:
        return self.kind == "uniform" or (self.a == 1 and self.b == 1)

    def key(self) -> tuple:
        return ("uniform",) if self.is_uniform else ("beta", self.a, self.b)


@dataclass(frozen=True)
class Model:
    """Uniform or scaled Beta immigration; c_i is TOTAL immigration intensity.

    Defaults are illustrative dimensionless parameters, not a clinical fit.
    Retain x > ell. Time and length units must be consistent with c_f*x**alpha.
    """

    L: float = 1.0
    ell: float = 0.01
    c_f: float = 2.0
    alpha: float = 0.0
    a: float = 1.0
    b: float = 1.0
    c_i: float = 20.0
    exit: Exit = field(default_factory=Exit)
    immigration: Immigration = field(default_factory=Immigration)

    def __post_init__(self) -> None:
        for name in ("L", "ell", "c_f", "a", "b"):
            finite(getattr(self, name), name, positive=True)
        finite(self.alpha, "alpha")
        if finite(self.c_i, "c_i") < 0 or self.ell >= self.L:
            raise ValueError("require c_i >= 0 and 0 < ell < L")
        if not isinstance(self.immigration, Immigration):
            raise TypeError("immigration must be an Immigration specification")
        if not isinstance(self.exit, Exit):
            raise TypeError("exit must be an Exit specification")
        if self.exit.kind == "pfb" and self.exit.boundary > self.L:
            raise ValueError("PFB boundary must not exceed L")
        if not math.isfinite(self.a + self.b):
            raise ValueError("a+b must be representable")
        for x in (self.ell, self.L):
            r = self.rho(x)
            e = self.exit.rate(x, self.L)
            if not math.isfinite(r + e):
                raise ValueError("combined rate overflows; rescale the model")

    def rho(self, x: float) -> float:
        try:
            r = math.exp(math.log(self.c_f) + self.alpha * math.log(x))
        except OverflowError as exc:
            raise FloatingPointError("fragmentation rate overflows; rescale the model") from exc
        if not math.isfinite(r) or r <= 0:
            raise FloatingPointError("fragmentation rate under/overflows; rescale the model")
        return r

    def sample_root(self, rng) -> float:
        """Do not condition on the cutoff, clip endpoints or redraw a root."""
        if self.immigration.is_uniform:
            return float(rng.uniform(0, self.L))  # preserve the existing stream
        q = float(rng.beta(self.immigration.a, self.immigration.b))
        x = self.L * q
        if not (math.isfinite(q) and 0 < q < 1 and math.isfinite(x) and 0 < x < self.L):
            raise FloatingPointError("Beta immigrant root rounded to an endpoint; rescale or use less extreme shapes")
        return x

    def fragmentation_key(self) -> tuple:
        """Parameters that cannot change when reusing a realised skeleton."""
        key = (self.L, self.ell, self.c_f, self.alpha, self.a, self.b)
        return key if self.immigration.is_uniform else key + (self.immigration.key(),)
