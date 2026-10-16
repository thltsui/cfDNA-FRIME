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

## Measured first execution, 12 October 2026

Code 10ac6a33b9254afa82fcc58cdf1679aabe416677 passed 135 tests in 20.61 seconds.
The six analytic checks generated 1,800 independent populations; the largest
mean discrepancy was 2.48 Monte Carlo SE. The complete 4,000-proposal resource
check took 1.49 seconds and stopped before screening. Uniform-negative had 271
proposals above the 100,000-family cap; Beta-negative had 275. Neither positive
case exceeded it. The largest planned population required 31.64 million families,
and the full screen required an expected 763.59 million families.

These are prospective work counts, not generated populations or goodness-of-fit
results. Full records and validation uncertainty are in ignored
`results/immigration-sensitivity-20261012/`. A revised abundance prior has been
proposed to the user; the revised plan below was authorized on 14 October.

## Approved abundance prior, 14 October 2026

The full four-case screen uses mu=E|U_W| log-uniform on [100,10000] in every
case. W is the observable window above 39.5 bp, not the biological 1 bp window.
The intensity is total C_I=mu/K, where K is expected integrated observable
family count. The normalized immigration density includes roots below cutoff;
it is not conditioned on retention. This shared abundance prior replaces the
old rate prior for a new separately seeded plan; the old plan remains intact.

For uniform splitting let h=rho+E, V(x)=integral_ell^x v(y)dy. Then
v=1/h+2*rho*V/(x*h). Setting W=V/x^2 and t=log(x/ell) gives
dW/dt=1/(x*h)-2*E*W/h, W(ell)=0. Finally K=integral p_I(x)v(x)dx.
`abundance_prior.py` solves this stable scalar equation in log length, splitting
at the PFB boundary, and integrates the Beta density with endpoint-aware
quadrature. Every coefficient is checked with tenfold tighter tolerances;
quadrature estimates and refinement discrepancies are saved. They are numerical
diagnostics, not rigorous interval error bounds. SciPy is an optional `abc`
dependency. Analytic no-exit and constant-exit checks cover both root laws.

This calculation only defines the prior coordinate. Screening and confirmation
still generate independent random Poisson-forest populations with the existing
p-mass certificate. It is not a deterministic mean fit. The mean count is a
latent abundance/dispersion parameter, not an estimate of plasma concentration.
Engineering caps may be raised prospectively after the entire plan is inspected;
no row is filtered and any realized cap failure remains explicit.

## Completed full screen, 14 October 2026

Scientific source: `4d99ff47b2f8ae1b77fe3f7b38eb808fcefe74e7`. Frozen root seed
`2026101401`. All 4,000 settings passed the revised preflight; all 5,200 random
populations completed, with no resource failures, no-read outcomes, skipped rows
or replacements. Screening and confirmation streams are disjoint; all recorded
latent-population and read streams are unique. Source, plan, data, histogram and
summary consistency checks passed.

The shared prior is on expected observable count, not fixed total immigration
rate; this affects finite-population variability and is part of the scientific
specification. Its deterministic mapping used tolerance refinement and retained
all 4,000 draws. One difficult mapping required additional refinement; no draw
was filtered. The largest accepted relative refinement discrepancy was 9.03e-8.

The target has 5,572 reads with none excluded. The discrepancy D is the maximum
absolute difference of the read and observed CDFs at the fixed bin edges. The
table's displayed candidate is selected **after confirmation**, by the smallest
maximum gap of the confirmation mean CDF. This is a descriptive display rule,
not an independently validated model-selection procedure. All 12 candidates are
retained in the local tables.

| Case | Screen D<=.05/.10/.15, out of 1,000 | Display ID | alpha_F | Mean-CDF max gap | Median individual D | Fresh D<=.05 |
|---|---|---:|---:|---:|---:|---:|
| Positive, uniform | 0 / 2 / 17 | 115 | 0.01720 | 0.09451 | 0.11181 | 0/100 |
| Positive, Beta | 0 / 0 / 10 | 1303 | 0.08123 | 0.09635 | 0.12078 | 0/100 |
| Negative, uniform | 2 / 65 / 199 | 2532 | -0.63220 | 0.04000 | 0.04684 | 70/100 |
| Negative, Beta | 6 / 75 / 195 | 3269 | -0.96017 | 0.02437 | 0.03033 | 100/100 |

Marginal Wilson 95% intervals for the last column are [0,.0370] for 0/100,
[.6042,.7811] for 70/100, and [.9630,1] for 100/100. These intervals condition on
candidate parameters and do not account for the subsequent display selection.
Neither the threshold nor the screen acceptance ratio is a calibrated test or a
posterior model probability. With only six negative-Beta and two negative-uniform
screen acceptances at D<=.05, this is not a precise posterior inference for alpha.

Candidate 3269 has B=1,214.41 bp, e=0.92843, a_I=1.61794, b_I=2.95466,
mu=8,577.88 and total C_I=280.12894 in the rho(B)=1 gauge. Its mean predicted
tail fractions, observed fractions and Monte Carlo SEs are:

| Length | Observed | Predictive mean | MCSE of mean | Pointwise 95% predictive interval |
|---|---:|---:|---:|---:|
| >1 kb | .71877 | .73440 | .000804 | [.71855,.74881] |
| >5 kb | .24390 | .24026 | .000781 | [.22555,.25540] |
| >10 kb | .04864 | .05078 | .000354 | [.04494,.05746] |

There is still a systematic residual: at 2,520 bp, the mean CDF differs by
0.024370 with MCSE 0.000872. The observed CDF is below the pointwise 2.5%
predictive quantile there. Pointwise intervals are not simultaneous bands and
do not supply a calibrated goodness-of-fit p-value. The best negative-uniform
candidate has a heavier >10 kb tail (.08863).

Negative-Beta candidates 3243 (alpha=-1.34455) and 3034 (alpha=-1.48029) also
give close but imperfect predictions: mean-CDF gaps .03369/.04174 and 66/63 of
100 populations with D<=.05. Their conditional Wilson intervals are
[.5628,.7454]/[.5322,.7182]. The best displayed exponent is above -1; this study
does not establish alpha<-1 as a threshold or identify a biological mechanism.
The Beta input improves the best candidate, but screen rates 75/1000 versus
65/1000 at D<=.10 do not establish a decisive case-wide advantage.

### Execution and artifacts

The production stages ran from 08:12:06 to 11:22:38 UTC on 14 October:

| Stage | Wall seconds | Families | Generated nodes |
|---|---:|---:|---:|
| Screen | 7,301.40 | 138,423,518 | 13,619,278,301 |
| Confirmation | 4,129.57 | 85,092,678 | 7,797,688,713 |
| Total driver | 11,431.96 | 223,516,196 | 21,416,967,014 |

Eight persistent family workers, batch size 128, M4 Mac mini with 10 cores and
16 GiB RAM. Sampled aggregate peak RSS was 812,482,560 bytes (775 MiB), not an
exact OS high-water mark. Full suite: 158 passed in 21.69s; Python 3.14.3,
NumPy 2.5.3, SciPy 1.18.1. Per-population history delta=1e-10; the union bound
over this adaptively selected 5,200-population study is <=5.2e-7. This controls
omitted immigration history only, not observation error, ABC approximation,
prior-map numerical error or parameter uncertainty.

The complete plan was inspected before prospectively setting engineering caps
to 50 million families, 1 million nodes per family and 250 million generated
nodes per population. Ten resource-probe settings plus two 4/8-worker replays
were checked before production. Replays reproduced histograms and counts.
Resource caps were never reached by the production run.

All outputs are local and ignored under `results/immigration-sensitivity-20261014/`:
`abundance-plan.json`, `execution.json`, `tests.xml`, `validation-summary.json`,
resource-probe records, stage logs, and `full-screen/` with complete JSONL/NPZ
records. Analysis outputs are `full-screen/analysis/four-case-predictive-cdfs.png`,
`case-summary.csv`, `all-candidates.csv`, `all-candidate-tails.csv`,
`displayed-candidate-tails.csv`, `analysis.json` and `INTERPRETATION.txt`.
The local `analyze.py` recomputes these summaries without new simulations.

Reproduction uses the saved, fully specified plan and fresh output directories:

```sh
.venv/bin/python examples/immigration_sensitivity.py --stage check --plan results/immigration-sensitivity-20261014/abundance-plan.json --output results/reproduction
.venv/bin/python examples/immigration_sensitivity.py --stage screen --output results/reproduction --workers 8 --batch-size 128
.venv/bin/python examples/immigration_sensitivity.py --stage confirm --output results/reproduction --workers 8 --batch-size 128
```

Read [NEGATIVE_ALPHA_NOTE.md](NEGATIVE_ALPHA_NOTE.md) for the mathematical audit.
Negative alpha is allowed by the thesis's weak stationary framework. A singular
zero-cutoff density and dust flux require care; they do not invalidate these
positive-cutoff simulations. No production mathematical convention was changed.
