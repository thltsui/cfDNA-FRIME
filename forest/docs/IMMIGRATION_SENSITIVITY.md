# S04 immigration and fragmentation sensitivity

This authorized exploratory comparison asks whether flexible immigrant lengths
can explain the long-fragment profile without a negative fragmentation exponent.
It is not a reproduction of the thesis fits or a calibrated model-selection test.
The earlier pilot excluded all four thesis fitted model/parameter combinations.

## Frozen comparison

Four cases cross alpha in [-1.5,0) or [0,1.5] with uniform or Beta immigration.
In the Beta cases a_I and b_I are independently log-uniform on [0.5,8]. All cases
keep uniform splitting, PFB exit exponent -1, B log-uniform on [50,6000] bp,
and e=E(B/2) log-uniform on [0.5,4], in the time gauge rho(B)=1.
The original nuisance prior is total C_I log-uniform on [10,100]. Its negative
alpha corner can require tens of millions of historical families. The initial
preflight records this failure without filtering parameter draws. Any revised
prior requires an explicit scientific decision and a separate frozen plan.

The target is the largest profile, S04, at its full depth and length window.
L=16,569 bp and the exact observable projection removes lengths <=39.5 bp.
The observation law and fixed CDF bins match ABC_PILOT.md. Each positive-depth
readout is multinomial conditional on a random population. Empty populations
remain explicit no-read outcomes; they are never rerolled.

Screen 1,000 independent parameter draws per case, with one independent
population and read draw each. Select the three finite-distance settings with
smallest CDF gap per case, breaking ties by proposal ID. Then generate 100 fresh
populations and read draws per selected setting. This gives at most 5,200
populations. Screening and confirmation use disjoint seed namespaces.

Report full-CDF distances, fractions above 1/5/10 kb, and success frequencies at
the common descriptive distance thresholds 0.05, 0.10 and 0.15. Report Wilson
intervals across independent populations, predictive intervals and Monte Carlo
SEs. Best-of-1,000 screening distances are optimistic. Confirmation remains
in-sample and selected; case success ratios are not posterior model odds.

## Law, resources and provenance

Beta changes only the normalized root law, not total Poisson arrival intensity.
The p-mass certificate uses the conservative full Beta moment; see MATHEMATICS.md.
No adapted certificate or clock-integrated replacement for random populations is
introduced. Per-population delta is 1e-10; the planned 5,200-population union
budget is 5.2e-7, distinct from posterior approximation error.

Plan JSON contains all parameter draws and complete model specifications before
screening. Check records every lookback certificate and prospective family count.
Every simulation records code/data/plan provenance, model, streams, histogram,
read counts, certificate, generated nodes and runtime. Resource or numerical
failures stop the stage; partial outputs are diagnostic only. No costly or failed
proposal is silently replaced. Raw outputs remain in ignored local results/.

From forest/:

```sh
.venv/bin/python examples/plan_immigration_sensitivity.py --output results/sensitivity-plan.json
.venv/bin/python examples/immigration_sensitivity.py --stage check --plan results/sensitivity-plan.json --output results/sensitivity
# Only after the prospective resource check succeeds:
.venv/bin/python examples/immigration_sensitivity.py --stage screen --output results/sensitivity
.venv/bin/python examples/immigration_sensitivity.py --stage confirm --output results/sensitivity
```

The script refuses overwritten stages or a change in scientific source, input
data or frozen plan between stages. Worker count is an execution setting; the
family-indexed random streams do not depend on it.
