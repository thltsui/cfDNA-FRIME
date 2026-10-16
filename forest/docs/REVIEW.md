# FRIME LOCAL REPORT — 16 October 2026

Objective:
Confirm S01/S05 versus S06–S09 with fresh populations and Beta immigration, then assess uncertainty before testing a shared exponent.

What ran:
- 3,600 independent populations: six profiles, three prespecified candidates per exponent sign, 100 repeats each; finite destruction cutoff retained.
- 1,000 read bootstraps per profile and refits of all fresh readouts against the complete fixed Beta ABC bank.
- Provenance, stream independence, population counts and every saved histogram/distance audited.

Main result:
The candidate-fit contrast survives fresh population randomness. For the original rank-1 candidate of each sign, these are counts satisfying D<=.10:

| Profile | Positive candidate | Negative candidate |
|---|---:|---:|
| S01 | 91/100 | 48/100 |
| S05 | 55/100 | 16/100 |
| S06 | 12/100 | 89/100 |
| S07 | 28/100 | 66/100 |
| S08 | 0/100 | 99/100 |
| S09 | 0/100 | 100/100 |

These estimate q_theta=P_theta(D(Y_new,y)<=.10), not P(alpha_F<0|y). Read resampling also preserves the broad contrast, conditional on the existing bank.

Evidence:
- Runtime: 42.76 minutes; sampled aggregate peak RSS 458.6 MiB.
- Scale: 46,007,254 immigrant families; 4,692,970,673 generated nodes; eight workers, batch 128.
- Seeds/repeats: 2026101601; 100 fresh populations per candidate; no failed or replaced populations.
- Tests: 168 passed in 22.75s; count means differ from their analytic predictions by at most 3.04 empirical MCSE across 36 candidates.
- Commit: simulation source 32c07ff; new-study forward-history error bound <=3.6e-7.

Interpretation:
The contrast is reproducible among the tested parameter settings. Different biological exponents are not yet established. The identical alpha_F=-.814 candidate fits S06 in 100/100 and S08 in 99/100 draws, so their differing ABC medians need not mean different exponents.

Uncertainty / blocker:
Near-zero signs are poorly resolved: with true alpha_F=+.028 for an S06 candidate, the inferred median is negative in 88/100 repeats, while all 100 empirical ABC parameter ranges include the true value. Immigration/exit trade-offs, finite-bank resolution and tolerance remain important. No shared-exponent model has been tested.

Decision for me:
Proceed to a shared-alpha versus profile-specific-alpha comparison, allowing profile-specific Beta immigration and exit? I recommend a recovery-checked comparison, rather than treating the present medians as biological differences.

Details: [PROFILE_CONFIRMATION.md](PROFILE_CONFIRMATION.md). Local evidence: `forest/results/profile-confirmation-20261016/`. [Draft PR #1](https://github.com/thltsui/cfDNA-FRIME/pull/1) remains unmerged.
