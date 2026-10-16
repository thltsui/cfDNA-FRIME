# Independent population Monte Carlo; preparation for ABC

For a fixed parameter vector theta, the current engine samples

`U_(theta,T)^(1), ..., U_(theta,T)^(R) iid ~ Law(U_(theta,T))`.

Each replicate has independent Poisson immigration and independent immigrant
families. No realised tree bank is shared between replicates. Both lazy and
skeleton modes simulate random populations. This is distinct from a
clock-integrated estimator of the stationary mean measure.

For length bins B_j, retain H_rj = <1_(B_j), U_r>, abundance N_r=<1,U_r>, and
retained mass M_r=<x,U_r>. The interface returns these raw quantities for every
replicate, rather than retaining only their averages. For any recorded statistic S,

`mean(S) = sum_r S_r / R`,
`MCSE(mean(S)) = sample_sd(S_1,...,S_R) / sqrt(R)`.

The SE is estimated across independent whole populations. It is unavailable for
R=1. These estimates target the finite-history law. Sampling variation, omitted
history, and the eventual data observation model are separate sources of error.

## API and reproducibility

Use a script with a main guard on macOS:

```python
import numpy as np
from frime_forest import Model, Exit, MonteCarlo

if __name__ == '__main__':
    theta = Model(L=1, ell=.02, c_f=2, alpha=1, c_i=100, exit=Exit('pfb'))
    bins = np.linspace(theta.ell, theta.L, 51)
    with MonteCarlo(workers=4, batch_size=128) as mc:
        first = mc.run(theta, 1000, bins=bins, seed=20260914, delta=1e-9)
        more = mc.run(theta, 1000, bins=bins, seed=20260914,
                      start_replicate=1000, delta=1e-9)
    print(first.statistics())
```

The same `(seed, replicate_id)` always reconstructs the same population, within
the same model, delta, mode and numerical environment. `n`, worker count, batch
size and previous calls do not determine the random stream. Reusing IDs replays
populations; use nonoverlapping ID ranges or a new root seed for fresh Monte Carlo.
Same IDs across parameter values deliberately couple the runs and must not be
mistaken for independent replicates. No tree is reused within a population.

A versioned integer namespace plus Cantor pairing encodes the root seed and
replicate ID injectively. The derived integer seed is recorded, so
`sample(theta, seed=result.seeds[r], delta=..., mode=...)` reproduces a row.
The existing family/component SeedSequence streams remain unchanged. As with
all PRNG simulation, ideal independence is implemented numerically; it is not a
literal theorem about finite deterministic PRNG outputs. See the
[NumPy parallel RNG documentation](https://numpy.org/doc/stable/reference/random/parallel.html).

The context keeps a spawn pool alive across calls and dispatches bounded batches
of immigrant families. Workers return histograms, counts and ordered per-family
mass scalars. The parent uses a fixed mass reduction, so summary results are
invariant to batch partition. Bins may cover only a subwindow; total count and
mass still describe all retained fragments. No normalization is performed.

Existing family/population caps remain engineering failures. If any replicate
exceeds a cap, the call raises with its ID, closes the context, and returns no
successful prefix. Never discard such failures and treat successful seeds as iid.
Explicit mathematical zero exits are valid; underflow of a positive PFB exit rate
now raises. Nonrepresentable summary reductions also raise.

## History tolerance for an ensemble

`delta` is a PER-REPLICATE ideal-law tolerance. For common model and certificate,
independent replicate couplings give

`d_TV(Law(U_T)^R, Law(U)^R) <= 1 - exp(-R * lambda_bound) <= R * delta`.

The result exposes the first upper bound as `joint_history_bound`. Deterministic
summaries inherit it. To guarantee an ensemble tolerance eta for a planned total
R_total, explicitly choose `delta <= eta/R_total` and retain it across chunks.
The implementation never silently changes delta based on the size of a call.
A TV bound alone is not an absolute bias bound for unbounded count/mass statistics.
The existing p-mass certificate and float64 caveat remain; the adapted certificate
has not been implemented.

## Recorded experiments

From `forest/`:

```sh
.venv/bin/python examples/monte_carlo.py --n 1000 --workers 4 --delta 1e-9 --output results/mc-001
```

This writes `metadata.json` and `populations.npz` in a fresh directory. They record
source commit/status/hashes, machine and package versions, model, length bins,
root and derived seeds, IDs, certificate, per-replicate counts/masses/histograms,
immigrant/generated-node counts, timings, means and MCSEs. The first replicate that submits work incurs pool startup; subsequent
replicates reuse that pool. `--model-json`
accepts Model fields with a nested Exit object; `--bins-json` accepts explicit
bin edges. Defaults remain illustrative dimensionless parameters, not an mtDNA fit.

## Eventual mtDNA comparison and ABC

The simulator provides the latent population law. An inference pipeline also needs
an observation kernel Q_phi(dy | U), so that

`U ~ pi_theta`, `Y | U ~ Q_phi(. | U)`.

Given observed mtDNA data y_obs, summaries S and discrepancy d, a rejection-ABC
target would have the form

`pi_epsilon(theta,phi | y_obs) proportional to prior(theta,phi)`
`  * P_(theta,phi){d(S(Y), S(y_obs)) <= epsilon}`.

Neither Q_phi, the prior, the summaries/discrepancy nor epsilon has been fixed by
this change. A possible conditional sequencing-count model is multinomial at
fixed read depth with size-dependent detection weights, but adopting it is a
scientific choice, not an implementation default. The prior, model adequacy,
ABC threshold and history tolerance are different objects.

Before fitting, agree on units and the retained length window, sequencing-depth
handling, detection/selection effects, and which sample/cohort is compared. A
normalized profile removes direct abundance information. Scaling all physical
rates together preserves the stationary law, so stationary observations alone
cannot identify the absolute shared timescale. Posterior predictive comparisons
will be needed to assess whether the fitted model reproduces the data.

## Existing mtDNA reference data (read only)

`data_preprocessing/mtDNA_modeling_data.csv` contains 68,811 frequency-table rows,
including 7,240 chrM rows across 13 sample IDs, all labelled nanopore. Retained
chrM lengths in this file range from 46 to 16,569 bp. Columns are an unnamed index,
id, chr, length, freq, study, condition, detail and specimen. These are aggregate
file/schema checks, not a biological analysis or a statement about representativeness.

`script/data_analysis.ipynb` source cell 9 defines `prob` as each frequency divided
by that sample's total frequency across chromosomes. Cells 11 and 12 then select
chrM and create 1-bp and 10-bp groups using `np.digitize(..., right=True)`, with edges
from 0 to 16,500. Later comparison plots renormalize individual profiles to their
own totals. A chrM probability relative to all cfDNA and a conditional mtDNA
length distribution therefore have different denominators.

The current NumPy histogram convention is left-closed/right-open except the last
bin; the historical grouping is right-closed. Integer observed lengths on bin
boundaries require an explicit common convention. The file includes lengths above
16,500, so the historical terminal grouping must not be silently treated as an
ordinary complete bin. A future data adapter should rebin the raw frequency table
under agreed edges/window, preserving zero bins and reporting excluded counts.
The original CSV, notebooks and preprocessing script were not modified.
