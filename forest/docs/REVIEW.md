# Review packet 01 — Poisson-forest FRIME

**Date:** 8 September 2026. **Status:** prototype implemented and locally tested;
not merged, not benchmarked on the user's Mac mini.

## Result

The new engine works independently of the old next-event programme. It represents
FRIME as independent immigrant families, applies inherited killing, stitches those
families onto Poisson immigration, and reports a deterministic lookback selected
from an analytic history-error bound. Existing source, notebooks and data remain
untouched. The source idea is thesis Section 8.1.2; the explicit Lyapunov certificate
and validation identities are separately derived in MATHEMATICS.md.

## Evidence

**60 tests passed** in the recorded local run; see TEST_RESULTS.txt. Tests include
ancestral pruning, inclusive/exclusive event boundaries, cutoff and mass checks,
parallel reproducibility, resource failures, analytic stationary first moments,
and independent chronological-reference comparisons of count, mass and squared
count for nonuniform Beta splits and all three exit functions.

Local host: Linux-6.18.35-x86_64-with-glibc2.41, Python 3.13.5, NumPy 2.3.5.
The tests ran with actual spawned workers, not mocked parallelism.

### Cold-start timing: a deliberately modest PFB workload

Model: L=1, ell=.02, c_f=2, alpha=1, c_i=100, a=b=1,
E(x)=max(x^(-1)-.4^(-1),0). Requested delta=1e-6.
Lookback T=70.067299; selected p=4.
Three seeds per setting; medians include process startup, IPC, tree work and
assembly but exclude the final snapshot query. No legacy-runtime comparison.

| Mode | Workers | Median elapsed |
|---|---:|---:|
| lazy | 1 | 0.345 s |
| lazy | 2 | 1.142 s |
| skeleton | 1 | 1.347 s |
| skeleton | 2 | 1.482 s |

For seed 12345, both modes simulated 6,997
immigrant families. Lazy mode stored 36,842
nodes versus 337,192 for complete skeletons.
Outputs were bit-for-bit identical across worker counts within each mode.
They need not be identical across modes.

**Interpretation:** early pruning is valuable under strong exit. This workload is
too small to amortise cold process startup; these measurements do not demonstrate
parallel speedup. Full banks remain useful when their cost is amortised over many
exit comparisons. Warm-pool reuse and larger batches are next performance tests,
not claimed features. Mac mini throughput remains unmeasured.

### Example run

Illustrative CON model, L=1, ell=.01, c_f=2, alpha=0, c_i=20, c_e=1,
horizon=2, seed=20260908, workers=2: T=14.273656;
323 immigrant families; 3779 stored nodes;
110 fragments at time zero. The analytic history mismatch bound
is 9.999995e-07. A small bound does not certify the
entire numerical implementation or remove Monte Carlo error.

## Your 15-minute review

Read MATHEMATICS.md Sections 2, 4 and 5. The scientific convention is total-rate
c_i, positive cutoff ell, independent families, and size-based exit. The proposed
acceptance is of the construction and bound, not a clinical fit or an asymptotic
performance claim. Leave any mathematical objection as a PR comment.

**Recommended next task:** benchmark on the Mac mini using the included runner,
then choose between tightening T and optimising amortised bank generation. Do not
change biological assumptions merely to improve throughput.

## Handoff

Only additive files in forest/ are part of this change. The next session should
first read PR feedback and this packet, then perform one approved task. No scheduled
agent, automatic merge, or unattended research run has been activated.
