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
