# FRIME LOCAL REPORT — 14 October 2026

Objective:
Check whether negative alpha_F is consistently supported across mtDNA profiles, allowing Beta immigration and retaining the finite destruction cutoff.

What ran:
- Screened the other 12 profiles against the complete 4,000-population bank; retained S04 as a reference.
- Generated 48,000 new multinomial read draws at each profile's actual depth; no new latent populations were needed.
- Compared positive/negative alpha with uniform/Beta immigration at the three frozen tolerances .05, .10 and .15; audited every distance and acceptance summary.

Main result:
Negative alpha is favoured in several profiles, but not all. At D<=.10 with Beta immigration, negative/positive acceptances are S01: 5/25; S05: 1/19; S06: 37/1; S07: 46/1; S08: 29/0; S09: 41/0. Each sign had 1,000 proposals. Their accepted-alpha medians are +.27, +.85, -.71, -.51, -.85 and -.82, respectively. S04 remains 75/0 with median -1.00. The broad S01/S05 versus S06--S09 contrast persists with uniform immigration and at D<=.15. At .05 most profiles have too few acceptances for useful alpha inference.

Evidence:
- Runtime: marginal screen and input audit 4.87 seconds; this reuses the previous 121.69-minute population screen.
- Scale: 13 profiles, 52,000 profile-level outcomes from 4,000 unique latent populations; new profiles contain 50--894 reads versus S04's 5,572.
- Seeds/repeats: new read seed 2026101402; original population seed 2026101401. No new candidate confirmations or repeated latent simulations in this extension.
- Tests: 161 passed in 22.59 seconds. All bank hashes, new distances, acceptance counts and alpha summaries audited; original S04 results preserved.
- Commit: screen code 6df5d25; latent-bank scientific source 4d99ff4. Approximate joint forward-history bound 4e-7 applies to the shared bank, not the conditional posterior.

Interpretation:
The screen suggests profile variation, not a demonstrated universal alpha or proven biological heterogeneity. The thesis slope relation beta=-alpha_F-2 permits negative alpha: its reported mtDNA slope -1.56 corresponds to alpha_F=-.44. Beta immigration generally modifies that relation, so slope alone no longer identifies alpha_F.

Uncertainty / blocker:
Sparse accepted samples, tolerance dependence, different read depths and possible observation-model mismatch. Cross-profile Monte Carlo errors are correlated through the shared bank. Accepted-alpha quantiles are finite-tolerance ABC summaries, not calibrated confidence intervals. No joint shared-alpha model was fitted.

Decision for me:
Should the next bounded study independently confirm the contrasting S01/S05 and S06--S09 regimes before fitting a shared-alpha model? I recommend yes.

Details: [CROSS_PROFILE_ABC.md](CROSS_PROFILE_ABC.md). Plots/tables: `forest/results/cross-profile-20261014/analysis/`. [Draft PR #1](https://github.com/thltsui/cfDNA-FRIME/pull/1) remains unmerged.
