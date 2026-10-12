# FRIME LOCAL REPORT — 12 October 2026


Objective:
Test whether random stationary FRIME populations can reproduce the 13 observed mtDNA length profiles through ABC.

What ran:
- Separate profile fits: 256 calibration draws and 2,048 independent production proposals, using a shared prior-predictive bank.
- One synthetic recovery check and 100 fresh posterior-predictive populations per profile, including the synthetic dataset.
- Repeatability check: 100 fresh populations at each of the ten closest settings for S04, the largest sample; eight sequencing draws per population.

Main result:
The ABC pipeline works, but close mtDNA alignment is not established. No production proposal matched a real profile within a binned-CDF gap of 0.05. S04's best single match was 0.0817; repeating its ten closest settings left mean-CDF gaps of 0.092–0.129. In the synthetic check, the true boundary B=500 bp gave posterior median 517 bp and a 95% interval of 365–810 bp. The other three parameters remained weakly identified.

Evidence:
- Runtime: 15.27 minutes for inference, predictive checks and candidate repeats, using eight workers on the M4 Mac mini.
- Scale: 4,705 populations; 16,235,987 immigrant families; 1,702,018,055 generated nodes.
- Seeds/repeats: root 20261012 with disjoint stage/population streams; 100 fresh populations per predictive profile and per fixed candidate.
- Tests: 99 passed in 10.12 seconds; output consistency audit passed. Additional 600-population analytic projection check showed deviations within 2.18 MCSE.
- Simulation commit: 9398870c4a3e68fc15b2d4631792fb6a8c9994e1. Raw data unchanged. Per-population history delta=1e-10; whole-study coupling budget <=4.705e-7.

Interpretation:
The simulator and inference machinery are usable. The fit captures a length scale, but the fitted profiles tend to miss long fragments: in 10/13 samples, the observed fraction above 5 kb exceeds the estimated 97.5th predictive percentile. This is a descriptive check, not a formal model rejection. Broad predictive bands do not establish close agreement.

Uncertainty / blocker:
The pilot fixes PFB exit, uniform splitting/immigration and equal fragment detectability. ABC tolerances are coarse (0.128–0.241 for real profiles). The search was finite and restricted; its results do not rule out FRIME. Absolute biological timescales are not identified. The additional diagnostic had one startup failure before any population returned; the wrapper was fixed and the same seeds rerun, with both records preserved.

Decision for me:
Approve a bounded S04 sensitivity study allowing alpha_F<0 and weaker exit e<0.5? I recommend it because the most repeatable candidates lie near these prior limits, before introducing another biological mechanism.

Protocol: [ABC_PILOT.md](ABC_PILOT.md). Detailed local output: `forest/results/abc-pilot-20261012/`; plots: `predictive-cdfs.png`, `largest-profile.png`, `synthetic-recovery.png`. Earlier Monte Carlo and benchmark evidence remains in the corresponding ignored results folders and Git history.
