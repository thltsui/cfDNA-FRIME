# Daily research handoff — proposed, not activated

Suggested local review time: 20:30–20:45 Europe/London. A future scheduled runner
should produce the review packet beforehand. This document is not a scheduler,
and neither ChatGPT nor this repository is currently configured to run it daily.

## State

Current deliverable: validated Beta immigration and a frozen four-case S04
sensitivity runner; see REVIEW.md and IMMIGRATION_SENSITIVITY.md.
Scientific basis: MATHEMATICS.md. Code/test evidence: TEST_RESULTS.txt.
Completed: Mac mini baseline, iid Monte Carlo, restricted ABC pilot, Beta engine
validation and the complete sensitivity resource preflight. The latter stopped
before screening: 546/4000 planned settings exceed the family cap. No costly
proposals were filtered. Close mtDNA alignment remains unestablished.
The existing draft PR #1 includes completed work and stays unmerged.

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

## Current bounded task

The user authorized Beta-vs-uniform immigration crossed with both signs of alpha,
1,000 screening draws per case and fresh candidate repeats. Implementation and
validation are complete. The original rate prior fails the resource preflight.

A scientific decision is pending: whether to use a common expected observable
abundance prior log-Uniform(100,10000). This is a new prior, not a change of time
units. Do not silently adopt it, alter the parameter ranges, filter costly draws,
or increase work caps to launch the original sweep. If authorized, implement and
validate its deterministic prior map, then preflight the entire revised plan.
That mean calculation would parameterize immigration only; simulations must
still produce random FRIME populations. It does not guarantee a cheaper run.

Keep the current observation law explicitly provisional. The earlier pilot
excluded all four original thesis fitted model/parameter combinations and is not
a replication or refutation of them. Adapted history certificates remain
unapproved research notes. No recurring runner has been activated.
