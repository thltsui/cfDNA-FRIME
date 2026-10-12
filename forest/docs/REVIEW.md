# FRIME LOCAL REPORT — 12 October 2026

Objective:
Compare whether Beta immigration explains S04's length profile without a negative fragmentation exponent.

What ran:
- Implemented normalized Beta immigrant lengths, retaining total C_I and a conservative p-mass history certificate.
- Ran the complete suite, including 1,800 independent populations for six analytic mean checks and worker/batch reproducibility checks.
- Froze 1,000 proposals for each of four cases and calculated every history certificate before screening.
- Updated draft PR #1 with the latest completed work; no merge.

Main result:
Beta immigration is implemented and validated, but the four-case mtDNA screening has not run. Its original rate prior fails the prospective resource check: 271/1,000 negative-alpha uniform cases and 275/1,000 negative-alpha Beta cases exceed 100,000 historical families. The worst proposal requires about 31.6 million families; the complete screen would require an expected 763.6 million. The entire plan stopped before population sampling, with no proposals filtered. This is a computational limitation, not evidence of poor mtDNA fit.

Evidence:
- Runtime: full suite 20.61 seconds; all-proposal preflight 1.49 seconds on the 10-core M4 Mac mini, 16 GiB RAM.
- Scale: 1,800 analytic-validation populations, 153,247 families and 1,057,901 generated nodes; largest analytic discrepancy 2.48 Monte Carlo SE.
- Seeds/repeats: frozen screening root 2026101202; validation has six fixed independent 300-population ensembles. Detailed seeds, means and uncertainty are recorded locally.
- Tests: 135 passed; Python 3.14.3, NumPy 2.5.3.
- Code commit: 10ac6a33b9254afa82fcc58cdf1679aabe416677. Clinical data and legacy code unchanged.

Interpretation:
The old immigration prior is expressed per boundary-fragmentation time. With negative alpha, large fragments split much more slowly; maintaining the same arrival rate can make population histories extremely costly. Raising caps or discarding those settings would not be an acceptable bounded study.

Uncertainty / blocker:
The screening prior needs an explicit scientific decision. A common prior on expected observable abundance is a possible alternative, but changes finite-population variability and still needs its own complete resource preflight. No new posterior or biological conclusion is available. The earlier restricted ABC pilot did not reproduce the original thesis fitted regimes.

Decision for me:
Use a common expected observable population-size prior, log-Uniform(100,10,000), for all four cases? I recommend this interpretable comparison, subject to a fresh resource check before any production simulations.

Protocol: [IMMIGRATION_SENSITIVITY.md](IMMIGRATION_SENSITIVITY.md). Local records: `forest/results/immigration-sensitivity-20261012/`. Prior ABC results remain in `forest/results/abc-pilot-20261012/` and Git history.
