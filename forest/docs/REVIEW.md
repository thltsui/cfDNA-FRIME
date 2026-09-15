# FRIME LOCAL REPORT — independent population Monte Carlo

**Date:** 14 September 2026. **Objective:** provide iid FRIME population Monte Carlo as the simulator component of eventual ABC against mtDNA data.

## What changed

- Added `MonteCarlo` / `monte_carlo`: independent replicate streams, persistent family workers, raw per-population histograms/counts/masses, and replicate-level Monte Carlo SEs.
- Kept the existing finite-history law and p-mass certificate. The interface reports the ensemble history bound separately; changing chunk size does not change a replicate.
- Added an executable example saving JSON/NPZ provenance. Workers transfer summaries rather than complete genealogies.
- Fixed the previously observed positive-exit underflow to raise explicitly. Exact zero exits remain valid. No clinical or legacy source changed.

## Evidence

All **84 tests passed in 4.87 s**. Checks include direct full-population replay, worker/batch/order/chunk invariance, empty populations, explicit failure on caps, and analytic moments with L=2 to detect immigration-normalization mistakes.

For 20 PFB populations per run, median total time over three seeds was **3.495 s** with the original serial full-population interface, **1.071 s** with four persistent summary workers, and **0.743 s** with ten. Ranges: 3.476–3.495, 1.066–1.154, and .742–.882 s. Histograms, counts and masses match the recorded reference streams exactly. The first parallel ensemble includes pool startup; later ensembles reuse the pool. Comparison includes different output/IPC paths, not an isolated speedup from worker count alone.

Three independent 1,000-population ensembles in each mode were compared with the analytic uniform-split, constant-rate profile. Across bin/count/mass summaries, max absolute standardized error was **3.743 MCSE**, below the predeclared six-MCSE threshold. Modes share root seeds and are coupled; there are three independent ensembles per mode, not six mutually independent ensembles. This is validation evidence, not a proof of law equality.

Detailed results and reproduction driver: ignored `forest/results/monte-carlo-20260914/`. Source hashes/status are recorded because these measurements preceded the commit. The mathematical/API guide is [MONTE_CARLO.md](MONTE_CARLO.md).

## Interpretation and next scientific decision

For fixed theta and T, the API generates `U_1,...,U_R iid ~ Law(U_(theta,T))`. Monte Carlo SE concerns replicate variation; `1-exp(-R*lambda_bound)` concerns history omission. The existing mtDNA file has 13 chrM sample profiles. Alignment to those observations requires an explicit observation kernel `Y|U`, length window and binning convention; sequencing read count is not automatically biological population abundance.

**For the first mtDNA comparison, condition on each sample's mtDNA read depth and assess its fragment-length distribution? I recommend this as the initial observation target; adopting a specific sampling/detection kernel remains a scientific decision.** ABC itself, priors and discrepancy thresholds are not implemented. The adapted certificate remains deferred pending explicit approval.
