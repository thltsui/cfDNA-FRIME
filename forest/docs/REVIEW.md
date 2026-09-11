# Review packet 05 — adapted-bound structural check

**Date:** 11 September 2026. **Scope:** one mathematical task; no simulator code changed.

## What I did

PR #1 remains open/draft with no comments, review submissions or review threads. The previous decision on implementing the adapted Lyapunov certificate is therefore still unanswered, so I did not modify `certificate.py` or start performance optimisation. I instead checked a structural point that would matter directly to a safe implementation: whether the proposed `max{1, ...}` ODE can switch branches repeatedly.

## Main result

For uniform splitting and `alpha=1`, define

`D_kappa(x)=C_F x + E(x) - kappa`, `H(x)=integral_ell^x h(y)dy`, and `g(x)=2 C_F H(x)/D_kappa(x)`.

On the active ODE branch, `H'=g`, so

`d log(g)/dx = (C_F - E'(x))/D_kappa(x)`.

The supported CON, PNB and PFB exit rates are nonincreasing, hence `E'(x)<=0` wherever differentiable. Therefore, whenever `D_kappa>0`, `g` is strictly increasing after it first reaches 1. The candidate has **at most one branch switch**: `h=1` initially, then the ODE branch forever. PFB's boundary is only a continuous kink and does not change this conclusion.

For the benchmark PFB case, the feasibility ceiling is also analytic: `inf_x(2x+E(x))=2B=0.8`, explaining the earlier numerical condition `kappa<0.8`. The below-boundary ODE is rational and has an elementary closed form. Re-evaluating it gives `bar_h=356.53609374898` and `T=34.2022504` at `kappa=.72`, reproducing the previous numerical result without an ODE solver.

Full derivation: [`TAIL_BOUND_NOTE.md`](TAIL_BOUND_NOTE.md).

## Evidence / uncertainty

This is an analytic derivation plus high-precision arithmetic for the displayed decimals. No simulator source changed, so repository tests were not rerun. The general additive-weight certificate remains unimplemented and the float64-versus-formal-interval caveat remains.

## One decision

**Approve implementing the adapted certificate for uniform splitting with `alpha=1`? I recommend yes: the single-switch proof removes the main structural ambiguity I would want resolved before coding it.**
