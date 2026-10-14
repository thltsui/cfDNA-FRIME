# Daily research handoff — proposed, not activated

Suggested local review time: 20:30--20:45 Europe/London. No recurring runner has
been activated; this file is not a scheduler.

## State

Completed: baseline, iid Monte Carlo, restricted ABC pilot, Beta immigration,
the full S04 sensitivity study (4,000 screening + 1,200 fresh candidate
populations), negative-alpha/cutoff audit, and the authorized marginal ABC screen
of all other 12 mtDNA profiles. The latter reuses the complete 4,000-population
bank and draws 48,000 new readouts; it does not generate new latent populations.

Current report: REVIEW.md. Cross-profile protocol/results: CROSS_PROFILE_ABC.md.
S04 protocol/results: IMMIGRATION_SENSITIVITY.md. Production law: MATHEMATICS.md.
Cutoff/asymptotic audit: NEGATIVE_ALPHA_NOTE.md. Validation: TEST_RESULTS.txt.
Raw outputs are ignored under results/. Draft PR #1 remains open and unmerged.

## Per-session instruction

Read the current draft PR/comments and latest report before new work. Complete
at most one approved task. Preserve the model; propose scientific changes before
implementing them. Validate relevant changes and report actual evidence. Keep
one concise review packet, detailed logs separately, and no growing backlog.
Never merge automatically, alter legacy/data files, or claim biological causes
from descriptive ABC fits. Fail resource/numerical problems explicitly; never
replace failed draws or filter expensive realizations.

## Current bounded task

The requested other-profile ABC screen is complete. S01/S05 favour positive-alpha
candidates while S04 and S06--S09 favour negative candidates under the current
priors, both immigration laws and the .10/.15 tolerances. Several other profiles
are sparse or tolerance-sensitive. This is not a joint shared-alpha test or proof
of biological heterogeneity; shared bank reuse induces Monte Carlo dependence.

The user explicitly retained the finite destruction cutoff. The underlying
1 bp destruction assumption and exact 39.5 bp observation projection are kept
distinct. No integer-splitting or zero-cutoff model was introduced. The thesis
tail-slope relation permits negative alpha for uniform input; Beta immigration
adds a slope term and confounds that simple inference. Adapted certificates
remain unapproved.

Suggested next decision: independent candidate confirmation for contrasting
profiles before a shared-alpha/hierarchical model. Neither new confirmations
nor that model have been started. No recurring task has been activated.
