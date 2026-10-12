# Daily research handoff — proposed, not activated

Suggested local review time: 20:30–20:45 Europe/London. A future scheduled runner
should produce the review packet beforehand. This document is not a scheduler,
and neither ChatGPT nor this repository is currently configured to run it daily.

## State

Current deliverable: the completed exploratory mtDNA ABC pilot, described in REVIEW.md and ABC_PILOT.md; iid population Monte Carlo remains documented in MONTE_CARLO.md.
Scientific basis: MATHEMATICS.md, with thesis-derived assumptions distinguished
from subsequent derivations. Current code/test evidence: TEST_RESULTS.txt.
Completed: Mac mini baseline, Monte Carlo validation, and a bounded ABC pilot with synthetic recovery, fresh predictive checks and S04 candidate repeats. Close real-data alignment is not established. The next scientific decision is whether to expand the pilot parameter support for a bounded sensitivity study.

## Per-session instruction

Read the current FRIME draft PR, review comments, docs/REVIEW.md and
MATHEMATICS.md. Resolve human feedback before beginning new work. Complete at
most one approved task. Preserve the model; propose scientific changes before
implementing them. Run the relevant tests and report the actual environment.
Update one concise review packet with objective, change, evidence, uncertainty
and at most one decision for Terence. Keep raw logs separate. Work only on the
feature branch; never auto-merge, publish clinical claims, or alter legacy data.

If an unanswered scientific decision blocks progress, examine that question or
report the blocker. Do not silently choose an answer or generate another backlog.
Do not claim that tests, experiments or scheduling ran unless execution records
exist. Resource caps must cause explicit failure, not biased partial samples.

## Near-term queue (not a promised timetable)

1. Obtain a decision on a bounded S04 sensitivity study across alpha<0 and e<0.5.
2. Predeclare its ranges and certificate/work budgets before simulation; do not silently filter costly proposals.
3. Keep the current observation kernel explicitly provisional; interpret sensitivity separately from a new posterior or biological validation.
4. Implement the adapted history certificate only after explicit approval.
5. Keep clock-integrated mean estimation separate from random-population Monte Carlo.

The current bounded task is complete; the next scientific sensitivity study awaits a decision. Scheduling can be activated only
through a real runner with verified repository access and a confirmed schedule.
