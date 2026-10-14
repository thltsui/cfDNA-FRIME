# Daily research handoff — proposed, not activated

Suggested local review time: 20:30–20:45 Europe/London. This document is not a
scheduler; no recurring runner has been activated.

## State

Completed: Mac mini baseline, iid Monte Carlo, restricted ABC pilot, Beta
immigration validation and the full four-case S04 sensitivity study. The user
approved the common expected observable abundance prior on 14 October after the
original rate-prior plan stopped at resource preflight. The revised frozen plan
completed all 4,000 screening and 1,200 candidate populations, with no failed or
replaced proposals. Scientific source was fixed at 4d99ff4 throughout.

Read REVIEW.md for the concise report, IMMIGRATION_SENSITIVITY.md for protocol
and measured results, TEST_RESULTS.txt for validation, and NEGATIVE_ALPHA_NOTE.md
for the zero-cutoff mathematical audit. MATHEMATICS.md remains the production
model specification. The draft PR #1 is open and unmerged.

## Per-session instruction

Read the current draft PR and review comments before beginning new work. Complete
at most one approved task. Preserve the model; propose scientific changes before
implementing them. Run relevant tests and report the actual environment. Update
one concise review packet with evidence, uncertainty and at most one consequential
decision. Keep raw outputs in ignored results/. Work on the feature branch; never
auto-merge, publish clinical claims or alter legacy data.

Do not claim experiments or scheduling ran without execution records. Resource
caps must fail explicitly; never filter expensive or failed proposals and treat
the remainder as iid. Do not create a backlog of unreviewed implementations.

## Current bounded task

The approved full screen is complete. The best displayed candidate has negative
alpha and Beta immigration, but retains a systematic CDF residual. The study
does not establish a biological mechanism, posterior model odds, or a threshold
at alpha_F=-1. Display choices were selected using confirmation results; these
are descriptive in-sample comparisons.

Keep the observation law provisional and distinguish this study from the thesis
fits, whose exit regimes differ. Finite positive cutoff remains part of the
model. Zero-cutoff PDE asymptotics are an audited research note, not a change to
the simulator or its history certificate. Adapted certificates remain unapproved.
The next proposed scientific decision is validation across additional profiles;
no new study or recurring task has been started.
