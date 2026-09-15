# Daily research handoff — proposed, not activated

Suggested local review time: 20:30–20:45 Europe/London. A future scheduled runner
should produce the review packet beforehand. This document is not a scheduler,
and neither ChatGPT nor this repository is currently configured to run it daily.

## State

Current deliverable: iid population Monte Carlo with persistent family workers, described in REVIEW.md and MONTE_CARLO.md.
Scientific basis: MATHEMATICS.md, with thesis-derived assumptions distinguished
from subsequent derivations. Current code/test evidence: TEST_RESULTS.txt.
Completed: Mac mini baseline and Monte Carlo validation. Next scientific task: agree the mtDNA observation target and length/normalization conventions before ABC.

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

1. Agree the mtDNA observation model and summary/discrepancy target.
2. Validate the resulting data adapter against the unchanged frequency tables.
3. Design ABC priors and synthetic-data recovery checks before clinical fitting.
4. Implement the adapted history certificate only after explicit approval.
5. Keep clock-integrated mean estimation separate from random-population Monte Carlo.

The immediate task is implementation review. Scheduling can be activated only
through a real runner with verified repository access and a confirmed schedule.
