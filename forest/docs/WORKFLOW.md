# Daily research handoff — proposed, not activated

Suggested local review time: 20:30--20:45 Europe/London. No recurring runner has
been activated; this file is not a scheduler.

## State

Completed: baseline, iid Monte Carlo, restricted ABC pilot, Beta immigration,
the full S04 sensitivity study (4,000 screening + 1,200 fresh candidate
populations), negative-alpha/cutoff audit, and the authorized marginal ABC screen
of all other 12 mtDNA profiles. The latter reuses the complete 4,000-population
bank and draws 48,000 new readouts; it does not generate new latent populations.

Fresh six-profile confirmation is also complete; protocol: PROFILE_CONFIRMATION.md.

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

Completed the authorized fresh confirmation for S01/S05 versus S06--S09 with
Beta immigration: 3,600 independent populations, conditional read bootstraps,
known-truth recovery and full provenance audit. See PROFILE_CONFIRMATION.md and
REVIEW.md. Candidate-fit differences persist, but biological exponent differences
are not established; near-zero median signs are unreliable and a common candidate
fits both S06 and S08. No shared-alpha/hierarchical model has been started.

The finite destruction cutoff remains: underlying 1 bp destruction and exact
39.5 bp observation projection are distinct. No model change or adapted
certificate was introduced. A shared-versus-profile-specific exponent comparison
with nuisance flexibility and recovery checks is the next proposed decision.
No recurring task has been activated.
