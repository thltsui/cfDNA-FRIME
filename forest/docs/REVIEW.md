# Review packet 02 — certificate audit

**Date:** 8 September 2026. **Scope:** mathematical audit only; no simulator code changed.

## What I checked

PR #1 remains open/draft and had no comments or submitted reviews at the start of this session. I re-derived the stationary-history argument in `MATHEMATICS.md` §§4–5, including the Poisson thinning step, the path coupling, the p-mass Lyapunov drift, and the interval lower envelope used for `kappa_p`.

## Main result

I find the present certificate mathematically valid under its stated assumptions: `0 < ell < L`, finite rates on `[ell,L]`, independent immigrant families, and Chapter-7 total immigration rate `c_i`. In particular,

`Lambda_T = c_i E[(tau_X-T)_+]`

is the mean number of omitted old families still alive at time 0, so the common-driver mismatch probability is exactly `1-exp(-Lambda_T)`. If no omitted family is alive at 0, it cannot reappear, so the same coupling controls the whole future path; `T` truncates immigration history and never forces family extinction.

The implemented p-mass estimate is conservative but sound. For the benchmark PFB model (`L=1, ell=.02, c_f=2, alpha=1, c_i=100, B=.4, beta=-1`), the selected `p=4` gives `kappa≈0.47391664` and `T≈70.0673` for `delta=1e-6`, reproducing the recorded value.

## Evidence on conservatism

As an independent diagnostic, I simulated 4,000,000 single immigrant families directly from competing fragmentation/exit clocks for that PFB model. Median extinction age was about 1.94; the 99.9999% empirical quantile was about 21.99; no sampled family survived past 24.66. This is **not** a certified tail estimate, but it strongly suggests that the current `T≈70` bound leaves substantial optimisation headroom.

No repository tests were rerun because no code changed and the branch was not materialised in the execution container; the Monte Carlo diagnostic was standalone.

## One decision

**Should the next session prioritise deriving a materially sharper *certified* survival/lookback bound before Mac-mini performance optimisation, rather than treating the current conservative certificate as the baseline?**
