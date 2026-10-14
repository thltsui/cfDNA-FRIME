# FRIME LOCAL REPORT — 14 October 2026

Objective:
Test whether Beta immigrant lengths can explain S04's mtDNA profile without negative fragmentation exponents.

What ran:
- Four frozen cases: positive/negative alpha_F crossed with uniform/Beta immigration, 1,000 proposals each.
- Three candidates per case, each with 100 fresh random populations and read draws.
- Full tests, resource probes, provenance audit and a review of negative-alpha PDE convergence against the thesis.

Main result:
All 5,200 populations completed without failures or replaced proposals. The best displayed negative-alpha/Beta candidate has alpha_F=-0.9602 and a maximum discrepancy between its mean read CDF and the observed CDF of 2.44 percentage points. Corresponding best displayed discrepancies are 4.00 points for negative/uniform, 9.45 for positive/uniform and 9.64 for positive/Beta. All 100 fresh draws for the best candidate had individual CDF distance <=0.05. However, a systematic residual remains: at 2,520 bp the 2.44-point mean discrepancy has Monte Carlo SE 0.087 points. This is improved alignment, not a demonstrated adequate observation model.

Evidence:
- Runtime: 190.53 minutes total; sampled aggregate peak RSS 775 MiB, on the 10-core M4 Mac mini, 16 GiB RAM.
- Scale: 223,516,196 immigrant families; 21,416,967,014 generated nodes; eight persistent workers, batches of 128.
- Seeds/repeats: root 2026101401; 4,000 screening populations plus 12 x 100 fresh candidate populations. No streams reused.
- Tests: 158 passed in 21.69 seconds; Python 3.14.3, NumPy 2.5.3, SciPy 1.18.1.
- Scientific commit: 4d99ff47b2f8ae1b77fe3f7b38eb808fcefe74e7. Per-population history tolerance 1e-10; study union bound <=5.2e-7.

Interpretation:
Within the tested priors and observation law, negative-alpha cases found substantially closer candidates. Beta immigration improved the best confirmed candidate but did not remove the preference seen in this screen for negative alpha. Candidates below -1 also fit reasonably; the best is above -1. These results do not estimate a reliable posterior probability of alpha_F<-1 or identify a biological mechanism.

Uncertainty / blocker:
One target, restricted exit/splitting laws, selected candidates and in-sample comparisons. Display choices use confirmation results, so intervals are descriptive and do not adjust for selection. A singular zero-cutoff density need not invalidate equilibrium; dust flux and the topology of convergence require separate care. The implemented positive-cutoff model is unchanged. This study does not replicate or refute the original thesis fits.

Decision for me:
Should the next study check these candidate regimes against additional mtDNA profiles before narrowing alpha_F? I recommend yes, with the observation law and discrepancy fixed beforehand.

Details: [protocol/results](IMMIGRATION_SENSITIVITY.md), [negative-alpha audit](NEGATIVE_ALPHA_NOTE.md). Raw records and plots: `forest/results/immigration-sensitivity-20261014/`. [Draft PR #1](https://github.com/thltsui/cfDNA-FRIME/pull/1) remains unmerged.
