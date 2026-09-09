# Review packet 03 — sharper lookback candidate

**Date:** 9 September 2026. **Scope:** one mathematical task; no simulator code changed.

## What I did

PR #1 is still open/draft and currently has no comments or review threads. Because the previous packet's decision on whether to tighten the certificate is unanswered, I investigated that question rather than moving on to Mac-mini optimisation.

I derived a more general additive Lyapunov certificate. Instead of the current fixed power weight `h(x)=(x/ell)^p`, choose any `h>=1` whose one-family generator satisfies `G V_h <= -kappa V_h`. Then

`P(tau_x>t) <= h(x) exp(-kappa t)`

and, for uniform immigration,

`Lambda_T <= C_I * mean(h(X)) * exp(-kappa T) / kappa`.

This preserves the existing Poisson coupling exactly; `T` still truncates immigration history, never family lifetime.

## Main result

For the benchmark PFB model (`L=1, ell=.02, C_F=2, alpha=1, C_I=100, a=b=1, B=.4, beta=-1`), uniform splitting makes the optimal-weight condition a one-dimensional ODE. Using the deliberately fixed `kappa=0.72`, the resulting analytic weight has `mean(h(X)) ≈ 356.536` and gives

**`T ≈ 34.20` for `delta=1e-6`, versus the current `T ≈ 70.07`.**

So a better certificate can cut expected pre-zero immigrant families from about 7,007 to 3,420 before touching parallelisation. The derivation is in [`TAIL_BOUND_NOTE.md`](TAIL_BOUND_NOTE.md).

## Evidence / uncertainty

The Lyapunov argument is analytic; the benchmark decimal evaluation used high-precision numerical quadrature/algebra, not Monte Carlo. I did not change `certificate.py` and did not rerun repository tests. The quoted decimals are not formal interval arithmetic, matching the numerical caveat already attached to the current certificate.

This cheap ODE form currently applies to uniform splitting with `alpha=1`; the general additive-weight theorem is broader, but non-uniform Beta splitting would require a Volterra integral bound.

## One decision

**Approve implementing this adapted Lyapunov certificate for the uniform-split, `alpha=1` case before Mac-mini benchmarking? I recommend yes, because it approximately halves the dominant history window without changing the FRIME law.**
