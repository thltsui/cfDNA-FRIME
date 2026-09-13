# Review packet 07 — conservative adapted-certificate algorithm

**Date:** 13 September 2026. **Scope:** one implementation-design task; no simulator source changed.

## What I did

PR #1 remains open/draft with no comments or submitted reviews, so approval to implement the adapted Lyapunov certificate is still unresolved. I did not change `certificate.py` or start performance optimisation. Instead I addressed the remaining numerical concern from the previous review: how to implement the tighter certificate without trusting root finding or a numerical ODE near `D_kappa=0`.

## Main result

For uniform splitting with `alpha=1`, the existing geometric size partition gives a conservative **interval supersolution**. On each interval `[l,u]`,

`D_kappa(x)=C_F x + E(x) - kappa >= d = C_F l + E(u) - kappa`.

Therefore the exact adapted equation

`H'=max(1, 2 C_F H / D_kappa(x))`

is dominated by the constant-coefficient equation

`Hhat'=max(1, (2 C_F/d) Hhat)`,

whose interval update is available in closed form. Propagating these updates gives `H_*(L) <= Hhat(L)`, so `Hhat(L)/L` is a conservative prefactor for the same Poisson-history error bound. No branch-switch root, ODE solver or quadrature is required. The current p-mass method remains a fallback.

For the existing PFB benchmark at `kappa=.72`, `delta=1e-6`, 256 geometric intervals give `T=34.8347`, versus the exact adapted value `34.2023` and the current p-mass value `70.0673`: about a **50.3% reduction** in certified lookback while retaining endpoint-envelope conservatism.

Full derivation and partition-refinement table: [`INTERVAL_ADAPTED_BOUND.md`](INTERVAL_ADAPTED_BOUND.md).

## Evidence / uncertainty

The recurrence and comparison proof are analytic; the displayed benchmark values were evaluated directly in Python. No simulator source changed, so repository tests were not rerun. Float64 remains non-formal arithmetic, exactly as for the current certificate.

## One decision

**Approve implementing this interval-supersolution certificate for uniform splitting with `alpha=1`, keeping the current p-mass certificate as fallback? I recommend yes; this removes the main numerical-implementation concern without changing the FRIME law.**
