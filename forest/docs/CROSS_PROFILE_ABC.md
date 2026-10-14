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
