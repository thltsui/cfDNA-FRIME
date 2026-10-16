# Conservative interval implementation of the adapted Lyapunov bound

**Date:** 13 September 2026  
**Status:** implementation design only; no simulator source changed.

## Goal

The adapted certificate in `TAIL_BOUND_NOTE.md` is mathematically attractive but an implementation based on root finding plus numerical ODE/quadrature creates an avoidable numerical-certification problem near

$$D_\kappa(x)=C_Fx+E(x)-\kappa=0.$$

For uniform splitting and `alpha=1`, the same geometric partition already used by `certificate.py` gives a conservative recursion that needs neither root finding nor a numerical ODE solver.

## Interval supersolution

The exact minimal adapted cumulative weight satisfies

$$H_*'(x)=F(x,H_*)=\max\left\{1,\frac{2C_FH_*}{D_\kappa(x)}\right\},\qquad H_*(\ell)=0.$$

On a partition interval `[l_j,u_j]`, the supported exit functions are nonincreasing, so

$$D_\kappa(x)\ge d_j:=C_Fl_j+E(u_j)-\kappa.$$

Reject a candidate `kappa` if any `d_j<=0`. Otherwise define

$$a_j=\frac{2C_F}{d_j},\qquad \widehat F_j(H)=\max\{1,a_jH\}.$$

Then `F(x,H)<=Fhat_j(H)` for every x in the interval. Starting from any upper bound `Hhat(l_j)>=H_*(l_j)`, scalar comparison shows that the solution of

$$\widehat H'=\widehat F_j(\widehat H)$$

remains an upper bound throughout that interval. The update is analytic. Writing `Delta=u_j-l_j`:

- if `a_j H0 >= 1`, then `H1=H0 exp(a_j Delta)`;
- if `a_j H0 < 1` and `1/a_j-H0 >= Delta`, then `H1=H0+Delta`;
- otherwise, with `t*=1/a_j-H0`, `H1=(1/a_j) exp(a_j(Delta-t*))`.

Propagating this over the partition gives

$$H_*(L)\le \widehat H(L),\qquad
\bar h_*\le \widehat H(L)/L.$$

Hence

$$\Lambda_T\le \frac{C_I\widehat H(L)}{L\kappa}e^{-\kappa T}$$

is a conservative past-history bound for the same Poisson-forest law. `T` still truncates immigration history, not family lifetime.

## Why this is preferable for v0.1

This reuses the existing deterministic geometric partition and endpoint lower-envelope logic. It removes the need to locate the branch switch, solve the adapted ODE numerically, or integrate `h` separately. Refining the partition tightens the supersolution. Float64 is still not formal interval arithmetic, so the existing numerical caveat remains; outward rounding could be added later.

A safe implementation can also retain the current p-mass certificate and return the shorter of the two lookbacks, so unsupported parameter regimes are unchanged.

## Benchmark check

For the existing PFB benchmark

`L=1, ell=.02, C_F=2, alpha=1, C_I=100, a=b=1, B=.4, beta=-1, c_e=1`, with `kappa=.72` and `delta=1e-6`:

| geometric intervals | upper `H(L)` | sufficient `T` |
|---:|---:|---:|
| 128 | 1021.4105 | 35.6641 |
| 256 | 562.1842 | 34.8347 |
| 512 | 441.7359 | 34.4999 |
| 1024 | 395.6976 | 34.3470 |

The exact adapted calculation reported in `TAIL_BOUND_NOTE.md` is `T≈34.2023`, while the current p-mass certificate is `T≈70.0673`. Thus even the conservative 256-interval recursion cuts the certified lookback by about 50.3% without relying on an ODE solver.

These decimals were evaluated directly from the analytic interval recursion in Python. They are numerical evidence for the implementation design, not a replacement for source-level tests once the method is approved and coded.
