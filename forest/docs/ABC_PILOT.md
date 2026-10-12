# Exploratory mtDNA rejection ABC

This is a bounded pilot, not an approved biological observation model or a claim
that FRIME explains the data. It simulates random FRIME populations with the
existing lazy engine. It does not substitute the stationary mean measure.

## Target and observation model

Fit each of the 13 chrM profiles in `data_preprocessing/mtDNA_modeling_data.csv`
separately. Use the raw weighted lengths, not normalized legacy notebook output.
Let L=16,569 bp, and analyze rounded integer lengths >=40 bp. The continuous
window is W=(39.5,L]. All current chrM reads fall inside this window. Bin edges
are half-integers: log-spaced thresholds plus 100,200,500,1000,2000,5000,10000 bp.

For a stationary counting measure U_theta, write h_j=U_theta(B_j). Given the
observed window depth n, the working read model is

    Y | U_theta,n ~ Multinomial(n, h/sum(h)).

This assumes equal detectability per fragment inside W and sampling with
replacement. It is not mass-weighted detection. If sum(h)=0, return an explicit
no-read outcome and assign infinite distance from a positive-depth observation;
never redraw the population until nonempty. Thus n is a planned read depth with
a failure category, rather than conditioning the latent law on successful reads.
Length-dependent sequencing bias, duplicate reads, extraction effects and
between-sample biological structure are not established by the input table.

A biological cutoff of 1 bp is used only to describe the underlying model.
The computation removes fragments <=39.5 bp exactly: since descendants are
smaller than their parents and families do not interact, an excluded subtree
can never return to W. Under common immigration and potential node marks,
P_W U_ell = U_39.5 pathwise. Uniform immigrant roots are still drawn on (0,L)
at the original TOTAL intensity; roots below the projection cutoff are empty.
The certificate and logged abundance/mass concern this observable projection,
not the unobserved full population. Data in W cannot identify a lower cutoff.

## Explicit exploratory prior

Normalize lengths by L, set b=B/L, and fix the time gauge rho(b)=1:

    rho(x) = (x/b)^alpha,
    E(x) = e * max(b/x - 1, 0),
    C_F = b^(-alpha), C_E = e*b, C_I = nu.

All prior coordinates are independent:

| Coordinate | Prior |
|---|---|
| alpha | Uniform(0,2) |
| B in bp | log-Uniform(50,6000) |
| e=E(b/2) | log-Uniform(0.5,4) |
| nu, total immigration per boundary-fragmentation time | log-Uniform(10,100) |

Splits are Uniform(0,1), exit is PFB with exponent -1, and immigrants are
Uniform(0,L). These are consequential pilot restrictions. All rates cannot be
identified absolutely from a stationary profile; the gauge fixes that symmetry.
Finite-population compositional variability can still depend on nu.

## Rejection law and checks

Distance is max_j |F_Y(j)-F_observed(j)| on the fixed ordered bins. Separate prior
calibration draws set epsilon to the ceil(0.05*n_calibration)th ordered distance
for each dataset. Freeze it before generating an independent production bank.
Only production draws with distance <=epsilon enter the approximate posterior:

    pi_epsilon(theta|y) proportional to
        pi(theta) P_theta{d(S(Y),S(y)) <= epsilon}.

The calibration quantile controls acceptance frequency, NOT goodness of fit.
Report its absolute value and fresh predictive discrepancies. This is finite
summary/tolerance ABC, not the exact posterior for all raw fragment lengths.
A shared iid (theta,U) bank is reused across *separate marginal fits*; each gets
independent multinomial read randomness. Each marginal rejection law is valid,
but cross-profile Monte Carlo estimates are correlated. Never multiply masks
or interpret this as a common-parameter or hierarchical joint posterior.

A separate synthetic dataset uses (alpha,B,e,nu)=(1,500,1.5,30), depth1000 and
fresh streams. Truth is not used in calibration/acceptance. Report truth versus
marginal posterior intervals; one recovery diagnostic is not simulation-based
calibration or evidence of nominal coverage.

For each fitted profile, uniformly sample accepted parameter draws, then generate
fresh U and Y. Store all outcomes, including no-read failures. Descriptive PPCs
show the mean predictive CDF, pointwise intervals, and a simultaneous band based
on the 95th percentile of sup-distance from that center. This band and its tail
rank are in-sample diagnostics, conditional on nonempty read outcomes, not a
calibrated goodness-of-fit test or external validation. Report no-read frequency
separately. Posterior acceptance simulations are never reused as PPC draws.

## Error control and execution

Default per-population history tolerance is 1e-10, using the unchanged p-mass
certificate. Sum the actual Lambda bounds to bound any history mismatch across
the finite experiment. This is not a posterior TV bound: rare-event conditioning
can amplify approximation. Float64 and pseudorandom error are outside the bound.
Adapted certificates remain unimplemented and unapproved.

The pilot uses 100,000 maximum historical families, 1,000,000 nodes per family,
and an explicit 10,000,000 generated-node population work cap. Small B permits
large trees above the no-exit boundary, so the package's default two-million cap
is insufficient for this prior. Workers return summaries in batches; all cap or
numerical errors abort the stage. Partial logs are diagnostic only, never a
filtered scientific sample. Failed seeds are not skipped or retried.

From forest/ (use a fresh output directory):

    .venv/bin/python examples/abc_pilot.py --stage probe --workers 4 --output results/abc-probe
    .venv/bin/python examples/abc_pilot.py --stage infer --workers 8 --calibration 256 --proposals 2048 --output results/abc-pilot
    .venv/bin/python examples/abc_pilot.py --stage ppc --workers 8 --ppc 100 --output results/abc-pilot

The seed and scientific settings must match across inference and PPC. Code hashes
and bins must also match. JSON manifests record hardware, versions, input hash,
source commit/status/hashes and assumptions; JSONL records every simulated
population's parameters, streams, certificate, histogram, nodes and runtime.
NPZ files contain observations, proposal banks, distances and predictive draws.
Results are ignored by git; source data are never modified.

For ABC tolerance interpretation see Wilkinson (2013),
https://arxiv.org/abs/0811.3355. Finite tolerance is an explicit inference
approximation separate from the forest history approximation.
