# Fresh confirmation of contrasting mtDNA profiles

Authorized 16 October 2026: independently confirm S01/S05 versus S06--S09,
allowing Beta immigration, then assess sampling uncertainty before considering
a shared fragmentation exponent. No common-alpha/hierarchical fit is included.

## Frozen experiment

Keep the finite destruction cutoff, uniform splitting, PFB exit exponent -1,
normalized Beta immigration and all parameter/abundance conventions from
CROSS_PROFILE_ABC.md. Use exactly the three lowest screening-distance candidates
within each sign of alpha for each of the six profiles, with proposal ID breaking
ties. Recompute this rule from the complete saved screen and verify it matches
the saved selections. This freezes 36 profile/candidate cells before any new
population is simulated. Rank 1 per sign/profile is the primary comparison;
ranks 2/3 are prespecified robustness checks, not post-confirmation replacements.

Generate 100 independent populations per cell, then multinomial readouts at that
profile's observed depth: 3,600 fresh populations in total. Root seed 2026101601;
latent namespace replicate_seed(root,900); replicate index cell*100+repeat;
read namespace (root,index,901). Even repeated candidate IDs across profiles
receive independent populations. All new streams are disjoint from earlier
screening and S04 confirmation. Empty populations remain no-read failures.

Use eight persistent family workers, batches of 128, and unchanged caps:
50 million historical families, 1 million nodes per family, 250 million generated
nodes per population. Certify every cell and preflight the complete plan.
Resource/numerical failures abort explicitly; no failed seed or expensive setting
is replaced. History tolerance is 1e-10 per population, so the new-study union
bound is <=3.6e-7; including the 4,000-population selection bank gives <=7.6e-7.
These are forward history bounds, not ABC posterior error bounds.

## Analyses specified before confirmation

1. Report every candidate's distance distribution, mean CDF discrepancy,
   fractions D<=.05/.10/.15, marginal Wilson intervals, predicted tails and
   pointwise prediction bands. Primary contrast uses screen-rank-1 candidates;
   do not replace them with whichever looks best after confirmation. None of
   these selected, in-sample diagnostics is a calibrated goodness-of-fit test.

2. For CDF bin j, let P_j be the latent population proportion and F_j the
   multinomial read proportion at depth n. Conditional on nonempty populations,
   Var(F_j)=Var(P_j)+E[P_j(1-P_j)/n]. Estimate both terms from independent
   populations. This separates process variability from read sampling at fixed
   parameters; it does not include parameter uncertainty.

3. Refit every fresh readout as a pseudo-observed dataset using the **complete
   fixed 2,000-proposal Beta bank**, including both signs, at the same profile
   depth and tolerances .05/.10/.15. Report acceptance counts, inferred alpha
   quantiles, median-sign recovery and unavailable fits. Retain no-acceptance
   outcomes in denominators. This is a local inference diagnostic at selected
   candidate truths and a fixed Monte Carlo bank, not global simulation-based
   calibration or proof of nominal credible-interval coverage.

4. Generate 1,000 multinomial resamples of each observed histogram at its original
   depth, with read-bootstrap component 902 and the stable profile index. Refit
   the same complete Beta bank. Report sign/median sensitivity and missing fits.
   These resamples approximate read-sampling uncertainty around the observed
   empirical profile; they do not reproduce population-to-population biological
   variation. The fixed bank's Monte Carlo uncertainty is not removed by this.

5. Describe the bootstrap distribution of the difference between the mean alpha
   medians of S01/S05 and S06--S09. Count unavailable group contrasts explicitly.
   This is conditional read-resampling sensitivity, not a p-value or confidence
   interval for biological heterogeneity: profiles were selected for the observed
   contrast, nuisance parameters differ, and inference uses a finite shared bank.

The primary tolerance stays .10; .05/.15 remain visible. Negative-alpha
acceptance frequencies and recovered parameter distributions are conditional on
the existing prior, cutoff, exit/splitting and observation laws. Finite repeats
reduce simulation uncertainty but do not identify a causal biological exponent.

## Commands and records

From forest/:

```sh
.venv/bin/python examples/confirm_profile_groups.py --stage prepare \
  --screen results/cross-profile-20261014 --output results/profile-confirmation-20261016
.venv/bin/python examples/confirm_profile_groups.py --stage run \
  --output results/profile-confirmation-20261016 --workers 8 --batch-size 128
```

The run refuses overwritten stages or changes to source, plan, data or screen
inputs. Each population records model, certificate, streams, worker settings,
histogram/readout, families, nodes, runtime and timestamp. The frozen plan and
run-start record contain Git/source hashes and machine/Python/NumPy metadata.
Raw output stays in ignored results/. Clinical data and legacy files are untouched.

## What the confirmation estimates

For one frozen parameter vector theta and the fixed observed profile y, define

\[
q_\theta(\epsilon;y)=\Pr_\theta\{D(Y_{new},y)\le\epsilon\}.
\]

The confirmation success count is Binomial(100,q_theta), conditional on the
screened candidate and observed data. Its Wilson interval describes Monte Carlo
uncertainty about this **fit probability**. It is neither a confidence interval
for alpha nor the probability that the biological exponent has a given sign.
The ABC sign probability would instead require the prior integrals

\[
\Pr_\epsilon(\alpha_F<0\mid y)
=\frac{\int_{\alpha_F<0}q_\theta(\epsilon;y)\pi(\theta)\,d\theta}
{\int q_\theta(\epsilon;y)\pi(\theta)\,d\theta}.
\]

Repeating only selected candidates does not compute those integrals. Therefore
the new simulations are not pooled into the old rejection bank as if they were
additional iid prior draws. The full fixed Beta bank is used for the diagnostic
refits, with its finite resolution and Monte Carlo dependence made explicit.

Conditional on a nonempty population, the binned read CDF satisfies

\[
\operatorname{Var}(F_j\mid\theta)
=\operatorname{Var}(P_j\mid\theta)
+\mathbb E[P_j(1-P_j)/n\mid\theta].
\]

Fresh populations measure both terms. Resampling observed reads as
Y*~Multinomial(n,y/n) addresses perturbations around the observed empirical
profile only. A persistent bootstrap contrast does not by itself establish
different biological alpha values. Local pseudo-data recovery asks whether the
same finite-bank inference can recover known candidate signs and magnitudes;
poor recovery qualifies interpretation even when a fit-probability contrast is
statistically clear. No p-value for biological heterogeneity is computed here.

## Completed results — 16 October 2026

All 3,600 populations completed with no failures, replacements or no-read outcomes.
| Profile | Positive candidate | Negative candidate |
|---|---:|---:|
| S01 | 91/100 | 48/100 |
| S05 | 55/100 | 16/100 |
| S06 | 12/100 | 89/100 |
| S07 | 28/100 | 66/100 |
| S08 | 0/100 | 99/100 |
| S09 | 0/100 | 100/100 |

Primary rank was frozen from screening; all 36 candidates are retained below. Counts are out of 100. Brackets are marginal Wilson 95% Monte Carlo intervals for q_theta at .10.

| Profile | Sign | Rank | Proposal | Alpha | D<=.05 | D<=.10 [interval] | D<=.15 | Mean CDF gap |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| S01 | positive | 1 | 1864 | 0.1717 | 2 | 91 [0.838, 0.952] | 100 | 0.0627 |
| S01 | positive | 2 | 1761 | 0.2352 | 0 | 69 [0.594, 0.772] | 100 | 0.0617 |
| S01 | positive | 3 | 1222 | 1.0348 | 0 | 19 [0.125, 0.278] | 90 | 0.1158 |
| S01 | negative | 1 | 3778 | -0.1023 | 0 | 48 [0.385, 0.577] | 96 | 0.0822 |
| S01 | negative | 2 | 3544 | -0.0892 | 2 | 53 [0.433, 0.625] | 91 | 0.0814 |
| S01 | negative | 3 | 3538 | -0.0811 | 1 | 33 [0.246, 0.427] | 82 | 0.0955 |
| S05 | positive | 1 | 1638 | 1.3903 | 0 | 55 [0.452, 0.644] | 100 | 0.0908 |
| S05 | positive | 2 | 1893 | 0.6197 | 1 | 95 [0.888, 0.978] | 100 | 0.0662 |
| S05 | positive | 3 | 1595 | 0.3653 | 0 | 80 [0.711, 0.867] | 100 | 0.0784 |
| S05 | negative | 1 | 3544 | -0.0892 | 0 | 16 [0.101, 0.244] | 70 | 0.1006 |
| S05 | negative | 2 | 3038 | -0.3750 | 0 | 1 [0.002, 0.054] | 17 | 0.1705 |
| S05 | negative | 3 | 3156 | -0.2617 | 0 | 0 [0.000, 0.037] | 19 | 0.1652 |
| S06 | positive | 1 | 1594 | 0.0280 | 0 | 12 [0.070, 0.198] | 89 | 0.1156 |
| S06 | positive | 2 | 1324 | 0.1053 | 0 | 5 [0.022, 0.112] | 63 | 0.1096 |
| S06 | positive | 3 | 1042 | 0.0306 | 0 | 12 [0.070, 0.198] | 79 | 0.0880 |
| S06 | negative | 1 | 3252 | -0.9621 | 18 | 89 [0.814, 0.937] | 100 | 0.0452 |
| S06 | negative | 2 | 3967 | -0.8137 | 58 | 100 [0.963, 1.000] | 100 | 0.0330 |
| S06 | negative | 3 | 3023 | -0.4124 | 8 | 99 [0.946, 0.998] | 100 | 0.0460 |
| S07 | positive | 1 | 1594 | 0.0280 | 0 | 28 [0.201, 0.375] | 100 | 0.0997 |
| S07 | positive | 2 | 1324 | 0.1053 | 0 | 0 [0.000, 0.037] | 45 | 0.1458 |
| S07 | positive | 3 | 1345 | 0.1852 | 0 | 0 [0.000, 0.037] | 96 | 0.1211 |
| S07 | negative | 1 | 3738 | -0.1741 | 1 | 66 [0.563, 0.745] | 98 | 0.0830 |
| S07 | negative | 2 | 3817 | -0.4873 | 41 | 100 [0.963, 1.000] | 100 | 0.0302 |
| S07 | negative | 3 | 3980 | -0.3700 | 21 | 99 [0.946, 0.998] | 100 | 0.0529 |
| S08 | positive | 1 | 1324 | 0.1053 | 0 | 0 [0.000, 0.037] | 45 | 0.1325 |
| S08 | positive | 2 | 1248 | 0.0865 | 0 | 0 [0.000, 0.037] | 57 | 0.1390 |
| S08 | positive | 3 | 1792 | 0.0649 | 0 | 0 [0.000, 0.037] | 3 | 0.2239 |
| S08 | negative | 1 | 3967 | -0.8137 | 51 | 99 [0.946, 0.998] | 100 | 0.0265 |
| S08 | negative | 2 | 3001 | -0.8686 | 4 | 93 [0.863, 0.966] | 100 | 0.0533 |
| S08 | negative | 3 | 3195 | -0.6274 | 1 | 53 [0.433, 0.625] | 98 | 0.0705 |
| S09 | positive | 1 | 1792 | 0.0649 | 0 | 0 [0.000, 0.037] | 16 | 0.1826 |
| S09 | positive | 2 | 1945 | 0.2558 | 0 | 0 [0.000, 0.037] | 12 | 0.1581 |
| S09 | positive | 3 | 1736 | 0.2276 | 0 | 0 [0.000, 0.037] | 52 | 0.1439 |
| S09 | negative | 1 | 3717 | -0.8745 | 73 | 100 [0.963, 1.000] | 100 | 0.0347 |
| S09 | negative | 2 | 3019 | -0.9835 | 13 | 79 [0.700, 0.858] | 100 | 0.0370 |
| S09 | negative | 3 | 3601 | -0.9675 | 52 | 99 [0.946, 0.998] | 100 | 0.0319 |

### Interpretation and inference diagnostics

The primary candidate comparisons favour positive exponents for S01/S05 and
negative exponents for S06--S09. This is a repeatable fit contrast conditional
on selected models, not a posterior odds calculation or a biological heterogeneity
test. Rank-2/3 confirmations show that a single screening winner need not remain
best under replication; they do not replace the primary candidates.

At tolerance .10, all 1,000 read-bootstrap group contrasts are available and
positive. The 2.5/50/97.5% resampling quantiles of the difference between mean
ABC medians (S01/S05 minus S06--S09) are 1.074/1.239/1.425. These are conditional
sensitivity summaries, not biological confidence intervals. At .05 only 2/1,000
group contrasts are available; no precise contrast should be inferred there.

Known-truth recovery qualifies the interpretation. For the S06 positive candidate
with true alpha=.027975, only 12/100 recovered medians have the correct sign at
.10, but all 100 empirical 2.5--97.5% ABC ranges include truth. None assigns
>=.95 negative mass. Thus median-sign errors near zero do not imply confident
exclusion of truth. For the S05 candidate with true alpha=1.3903, every recovered
median is positive, but the median of recovered medians is .8624 and 82/100
parameter ranges contain truth. Signs and magnitudes have different reliability.
These are local checks at selected truths using a fixed bank, not global SBC.

Identical proposal 3967, alpha=-.813728, satisfies D<=.10 in 100/100 fresh S06
and 99/100 fresh S08 draws; all its nuisance parameters are identical. This is
compatibility evidence within the negative group, not a fitted shared-alpha model.
Different original accepted medians do not establish distinct exponents.

Full threshold-dependent recovery and the population/read variance decomposition
are saved in analysis/recovery.csv and analysis/variance.csv. The latter evaluates
the two variance components at each candidate's largest mean-CDF residual; it does
not include parameter uncertainty or uncertainty in the observed latent profile.
No new simulations were pooled into the original iid-prior rejection bank.

### Reproduction and audit

After confirmation completes, from forest/:

```sh
.venv/bin/python examples/analyze_profile_confirmation.py results/profile-confirmation-20261016 --stage bootstrap
.venv/bin/python examples/analyze_profile_confirmation.py results/profile-confirmation-20261016 --stage confirm
```

Stages refuse existing output directories. The completed read bootstrap ran while
the independent population stage was running; its source hash and inputs are saved.
The archived analysis script and committed example have identical scientific calculations.
The committed example moves two plot labels to avoid overlaps; plot-manifest.json
records the presentation change and regenerated image hashes.
Run metadata, per-population JSONL, histograms, readouts, source hashes, certificates,
analysis tables and figures remain ignored locally. The independent count audit
is reproducible with results/profile-confirmation-20261016/audit_counts.py.

Run wall time 2565.714171s; 46,007,254 families; 4,692,970,673 nodes; sampled aggregate peak RSS 480854016 bytes. Full suite: 168 passed in 22.75s. Maximum count-mean discrepancy 3.040955 empirical MCSE. New forward-history union bound <=3.6e-7; including the 4,000-population selection bank <=7.6e-7. Neither is an ABC posterior error bound.
