# Review packet 06 — adapted bound formally dominates p-mass

**Date:** 12 September 2026. **Scope:** one mathematical task; no simulator code changed.

## What I did

PR #1 remains open/draft with no comments or submitted reviews, so the requested approval to implement the adapted Lyapunov certificate is still unresolved. I therefore did not change `certificate.py` or begin optimisation. I checked a narrower risk instead: could the specialised adapted certificate ever give a worse lookback than the existing p-mass certificate in its supported regime?

## Main result

**No.** For uniform splitting and `alpha=1`, fix any admissible `kappa` and define

`F(x,H)=max{1, 2 C_F H / (C_F x + E(x) - kappa)}`.

The adapted cumulative weight solves `H_*'=F(x,H_*)`, `H_*(ell)=0`. Any other feasible weight `h_0>=1`, with `H_0'=h_0`, obeys `H_0'>=F(x,H_0)`. Because `F` is nondecreasing in `H`, scalar ODE comparison yields

`H_* <= H_0` and `h_* <= h_0` almost everywhere.

The current p-mass weight `h_p=(x/ell)^p` is feasible whenever

`kappa <= inf_x {(1-m_p) C_F x + E(x)}`,

because cutoff can only remove child contributions. Therefore, at every decay rate certified by a p-mass candidate, the adapted construction has no larger prefactor and hence no larger sufficient lookback `T`. Optimising `kappa` can only improve further.

This gives a clean implementation safety rule: in the target regime the new method can be an optional refinement with the current certificate retained as fallback. It changes neither the Poisson-family law nor the meaning of `T`.

Full proof: [`TAIL_BOUND_NOTE.md`](TAIL_BOUND_NOTE.md).

## Evidence / uncertainty

This is an analytic comparison argument; no Monte Carlo or new numerical claims are needed. No simulator source changed, so repository tests were not rerun. Numerical implementation still needs conservative handling near `D_kappa=0`; float64 is not formal interval arithmetic.

## One decision

**Approve implementing the adapted certificate for uniform splitting with `alpha=1`, while retaining the current p-mass certificate as fallback? I recommend yes: we now have both structural simplicity and a formal no-regression result.**
