# Cross-profile mtDNA ABC screen

Authorized 14 October 2026: retain finite destruction cutoff and examine whether
negative fragmentation exponents remain compatible across the other mtDNA
profiles, allowing Beta immigrant lengths. This is an exploratory comparison of
separate profile parameters, not a common-parameter or hierarchical fit.

## Frozen design

Reuse the **entire** 4,000-population prior bank from the completed S04 study:
1,000 settings in each positive/negative-alpha x Uniform/Beta-immigration case.
Do not reuse only S04-selected candidates or accepted populations. Source,
plan/data hashes, complete coverage, histogram/log consistency and distinct
population streams must pass before new draws. Original scientific source is
4d99ff47b2f8ae1b77fe3f7b38eb808fcefe74e7.

All original priors and the observation law remain fixed: alpha sign strata
[0,1.5] and [-1.5,0], Beta immigration shapes independently log-uniform [.5,8],
expected observable abundance log-uniform [100,10000], uniform splitting,
PFB exit exponent -1, B log-uniform [50,6000] bp, e log-uniform [.5,4], rho(B)=1.
Total immigration intensity remains C_I=mu/K. The biological destruction cutoff
is 1 bp; computation uses the exact projection above 39.5 bp. Neither is fitted.

For each of the other 12 profiles, generate one independent multinomial readout
per saved population at that profile's original read depth. New read root seed
is 2026101402; streams use the original proposal ID and component 800 plus the
stable zero-based profile index. Preserve the original S04 reads and distances
as a reference. This makes 48,000 **new read draws**, not 48,000 new FRIME
populations. Each marginal screen uses 4,000 independent prior-predictive
populations. Cross-profile results share latent populations and parameter draws.

For each profile and immigration law, report both sign strata separately at
epsilon=.05,.10,.15. The middle value .10 is fixed for the primary display;
retain all tolerances in machine-readable output. For descriptive alpha summaries,
pool equal sign strata, corresponding to a Uniform(-1.5,1.5) prior conditional
on that immigration law. Report accepted count, fraction negative, and empirical
2.5/50/97.5% alpha quantiles. These are finite-summary, finite-tolerance ABC
sample summaries, not calibrated confidence intervals. Flag fewer than 20
acceptances as sparse. A zero count does not prove zero posterior probability.

Record the top three screening candidates in every case/profile for possible
follow-up; this run does not label them independently confirmed. No new latent
population generation is required for the requested marginal screen. No failure,
empty population or expensive proposal may be filtered or replaced.

## What can be learned

Consistently more negative-case acceptances across profiles would support the
negative-alpha region **within these priors and observation assumptions**.
Different best-fit exponents alone do not establish biological heterogeneity.
Broad overlapping accepted distributions do not establish universality either.
Other profiles have only 50--894 reads, so sampling uncertainty differs strongly
from S04 (5,572 reads). Conditional abundance and immigration shape can confound
alpha; close discrepancy must be checked before interpreting its distribution.

Do not multiply acceptance masks or posterior estimates across profiles: this
would impose an unrequested common full parameter vector and misuse dependent
Monte Carlo outputs. A shared-alpha hypothesis with profile-specific immigration,
exit and abundance needs its own joint/hierarchical model and model checks.
This screen can motivate that study but cannot establish it.

The history approximation is inherited from the 4,000 unique latent populations
(sum of recorded bounds approximately 4e-7). Reusing them for readouts does not
multiply the number of omitted-history events. This is a forward-simulation bound,
not a TV bound for an ABC posterior after acceptance conditioning.

## Tail-slope interpretation

The thesis's Chapter 6 relation for uniform splitting and uniform immigration,
in the PFB region x>B where E(x)=0, is

\[
u(x)\propto x^{-\alpha_F-2},\qquad
\beta=\frac{d\log u}{d\log x}=-\alpha_F-2.
\]

It permits negative alpha: a declining density with slope -1.56 corresponds to
alpha=-0.44. The thesis explicitly reports a mean mtDNA slope -1.56 and negative
alpha around -.5 (FRIME_comp.tex, lines 404--409). Its gDNA discussion reports
the steeper mean slope -3.32, corresponding to positive alpha around 1. These
are historical results, not newly estimated slopes for the present screen.

For a general normalized immigrant density p_I, uniform splitting and E=0 above
B instead give the exact stationary mean density on B<x<L:

\[
u(x)=\frac{C_I}{C_F}x^{-\alpha_F-2}
\left[x^2p_I(x)+2\int_x^L yp_I(y)\,dy\right].
\]

Writing the bracket as Q(x), differentiation gives Q'(x)=x^2p_I'(x), hence

\[
\frac{d\log u}{d\log x}
=-\alpha_F-2+\frac{x^3p_I'(x)}{Q(x)}.
\]

Uniform input has p_I'=0 and recovers the thesis formula exactly. Beta input
generally does not. For Beta(a_I,b_I) scaled to L,
p_I'/p_I=(a_I-1)/x-(b_I-1)/(L-x), so the upper tail can steepen because of
immigration even with negative alpha. A slope alone then does not identify
alpha without accounting for immigration. There is no x->infinity regime with
bounded support L; one must specify the observed range and exit-free regime.

The formula concerns the stationary mean number density. A single population fluctuates, and under the read model E[h_j/N] generally differs from E[h_j]/E[N] (conditioning on nonempty populations where needed). Therefore a fitted slope in one finite read sample is not an exact deterministic measurement of alpha, even with uniform immigration.

The slope here is of **number density**, or counts in equal-width linear bins.
Counts in equal-width log bins acquire an extra factor x. A cumulative survival
curve is a different integral and, with finite L, does not generally have the
same power-law slope. These distinctions precede any biological interpretation.

## Execution

From forest/, use a new output directory:

```sh
.venv/bin/python examples/cross_profile_screen.py \
  --bank results/immigration-sensitivity-20261014/full-screen \
  --output results/cross-profile-20261014 \
  --seed 2026101402
```

The output freezes its plan before read generation and records code/environment,
input hashes, profile depths, source bank commit, stream scheme, all read arrays,
distances, acceptance summaries and selected candidates. New code concerns only
screening/analysis; the simulator and clinical data remain unchanged.

## Measured results, 14 October 2026

Screening implementation `6df5d25`, original latent-bank scientific source
`4d99ff4`; new read seed 2026101402. All 12 additional profiles were screened
successfully, with S04 preserved as reference. Input audit and read generation
took 4.8696 seconds. No new latent populations were generated: 48,000 new
multinomial read draws use the original 4,000 independent prior-predictive
populations, whose generation took 7,301.40 seconds. All 52,000 marginal outcomes
are retained. No outcomes or proposals were filtered/replaced.

At the predeclared primary display tolerance D<=.10:

| Profile | Reads | Uniform negative/positive accepts | Beta negative/positive accepts | Beta accepted-alpha median |
|---|---:|---:|---:|---:|
| S01 | 267 | 4 / 17 | 5 / 25 | +.274 |
| S02 | 236 | 12 / 2 | 17 / 3 | -.212 |
| S03 | 190 | 8 / 11 | 16 / 8 | -.052 |
| S04 | 5,572 | 65 / 2 | 75 / 0 | -1.000 |
| S05 | 577 | 1 / 13 | 1 / 19 | +.855 |
| S06 | 432 | 37 / 0 | 37 / 1 | -.712 |
| S07 | 836 | 48 / 0 | 46 / 1 | -.508 |
| S08 | 349 | 32 / 0 | 29 / 0 | -.850 |
| S09 | 894 | 38 / 0 | 41 / 0 | -.823 |
| S10 | 50 | 3 / 0 | 4 / 0 | -.205 |
| S11 | 142 | 20 / 0 | 16 / 0 | -.279 |
| S12 | 80 | 4 / 4 | 6 / 1 | -.239 |
| S13 | 220 | 1 / 0 | 0 / 0 | unavailable |

Each sign/immigration case has 1,000 prior proposals. Pooling signs only within
an immigration case gives the alpha summaries; immigration laws are not mixed
into model odds. Quantiles use NumPy's default linear sample-quantile convention;
with small accepted counts they are unstable and interpolation should not be
mistaken for precise posterior support boundaries. In particular, individual
opposite-sign accepted draws can lie beyond an interpolated 95% range.

At D<=.15, Beta negative/positive counts remain 31/78 for S01 and 10/62 for
S05, versus 120/6, 122/10, 117/4 and 133/3 for S06--S09. Thus their broad
contrast is not confined to a single tolerance. S02 becomes 63/53 and S03
58/57, illustrating sensitivity near alpha=0. At D<=.05, most profiles have
zero or only a few acceptances; a precise small-tolerance posterior is not
available from this bank. S13 has no Beta acceptance at .10 but 30 at .15;
this may reflect limited proposal coverage or model discrepancy, not proof of
an impossible model.

The result supports investigating profile variation, but does not establish a
shared exponent, exclude every common-alpha model, or prove biological
heterogeneity. Different read depth, latent abundance, immigration, exit and
observation effects all matter. These screens are not independent biological
replications merely because they have separate profile labels. No new fresh
candidate confirmation or hierarchical model has been run in this extension.

Validation: full suite 161 passed in 22.59 seconds on the M4 Mac mini. New tests
check the multinomial mean against an analytic uncertainty bound, no-read
failures, stream reproducibility/separation, correct acceptance handling and
rejection of incomplete or modified banks. An independent post-run analysis
recomputed all 52,000 CDF distances and all acceptance/alpha summaries; checked
read depths, bank/output hashes and preservation of S04; all passed.

Local ignored artifacts: `results/cross-profile-20261014/plan.json`,
`screen.npz`, `case-summary.json`, `alpha-summary.json`, `selected.json`,
`screen-complete.json`, `validation.json`, `REPORT.txt`; `analysis/` contains
`alpha-summary.csv`, `case-acceptance.csv`, `audit.json` and figures
`alpha-epsilon-0.05.png`, `alpha-epsilon-0.10.png`, `alpha-epsilon-0.15.png`.
The retained `analyze.py` regenerates the static figures and independently audits
the recorded arrays. The machine-readable plan contains complete provenance and
the source ID mapping; public-facing plots use only S01--S13 labels.
