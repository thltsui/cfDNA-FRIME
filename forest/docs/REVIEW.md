# Review packet 04 — adapted-bound robustness

**Date:** 10 September 2026. **Scope:** one mathematical task; no simulator code changed.

## What I did

PR #1 is still open/draft and has no comments or review threads. The previous decision on implementing the adapted Lyapunov certificate is therefore still unanswered, so I did not change `certificate.py` or proceed to Mac-mini optimisation. Instead I tested whether the proposed bound depends sensitively on the arbitrary choice `kappa=0.72`.

## Main result

For the benchmark PFB model (`L=1, ell=.02, C_F=2, alpha=1, C_I=100, a=b=1, B=.4, beta=-1`), I independently re-evaluated the ODE weight and numerically minimised its certified lookback over feasible `kappa`.

The optimum is approximately

`kappa*=0.7262776`, `bar_h=441.3875`, `T=34.18862`.

The original simple choice `kappa=0.72` gives `T=34.20225`, only about `0.014` time units longer. So the earlier ~51% improvement over the current p-mass certificate (`T≈70.07`) is robust; it is not an artefact of tuning kappa.

There is also a useful stability warning: this model requires `kappa<0.8`. As `kappa` approaches `0.8`, the ODE denominator becomes small near the PFB boundary and the prefactor explodes; e.g. `bar_h≈2,277` at `kappa=.76` and `≈15,660` at `.78`, worsening the final lookback despite the faster exponential rate.

Full calculations are in [`TAIL_BOUND_NOTE.md`](TAIL_BOUND_NOTE.md).

## Evidence / uncertainty

The check used independent root finding, numerical quadrature and scalar optimisation. No Monte Carlo was used. No simulator source changed, so repository tests were not rerun. Decimal values are ordinary floating-point numerical evidence, not formal interval arithmetic.

## One decision

**Approve implementing the adapted certificate for uniform splitting with `alpha=1`? I recommend yes: the gain is large, the chosen kappa need not be finely tuned, and the simulator law remains unchanged.**
