# FRIME Poisson Forest — experimental v0.1

A new simulator implementing the immigrant-family construction proposed in
Section 8.1.2 of Terence Tsui Ho Lung's DPhil thesis. It does not import or modify
the legacy simulator or notebooks. Start with [the review packet](docs/REVIEW.md)
and [the mathematical specification](docs/MATHEMATICS.md).

## Implemented

Independent binary Beta fragmentation families; inherited size-dependent exit
pruning; Poisson immigration; an explicit moment-based stationary history bound;
reproducible process-parallel generation; frozen time queries; and three exit
families (CON, PNB, PFB). Uniform immigration on [0,L] is the current scope.

**Two modes:** `lazy` generates only the genealogy needed up to the observation
horizon and prunes early. `skeleton` generates complete size-truncated pure trees
and then prunes them. Both give the same ideal finite-history law, not necessarily
the same realisation for the same seed. For repeated exit comparisons, use the
explicit `presample` / `assemble` interface.

**Accuracy:** `delta` bounds ideal-law past-history truncation via the supplied
certificate. It does not bound numerical error, cutoff/model error or Monte Carlo
standard errors. No empirical burn-in or forced maximum family age is used.

## Install and run

From the repository root, using Python 3.10 or later:

```sh
cd forest
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[test]'
python -m pytest -q
python examples/demo.py --workers 2 --mode lazy
python benchmarks/run.py --workers 1 2 4 --repeats 3
```

The demo writes a snapshot, time-indexed histograms, seeds and certificate to
`results/demo/`. Its default model is illustrative and dimensionless, **not a
clinical parameter fit**. Examples use a `__main__` guard so process spawning
works in a terminal on macOS. Start with `workers=1` in an interactive notebook.
Worker count is explicit rather than assuming a particular Mac mini chip.

## Direct API

```python
from frime_forest import Model, Exit, sample
model = Model(L=1, ell=.01, c_f=2, alpha=1, c_i=100, exit=Exit('pfb'))
forest = sample(model, delta=1e-6, horizon=2, seed=123, workers=1)
sizes_now = forest.at(0)
sizes_later = forest.at(2)
print(forest.plan.certificate)
```

`c_i` means TOTAL immigration intensity, following Chapter 7. Some Chapter 6
formulae use `c_i L`; convert deliberately when reproducing old settings.
CON uses `Exit('constant', c_e=...)`; PNB/PFB power exponent defaults to -1.
For the thesis PNB/PFB normalisation use c_e=1. `ell` is a length cutoff, not a
KS convergence tolerance. No clinically calibrated units are supplied.

## Pre-sample and stitch

In a script, inside its `if __name__ == '__main__':` block:

```python
from frime_forest import make_plan, presample, assemble
plan = make_plan(model, delta=1e-6, horizon=2, seed=123)
bank = presample(model, len(plan.arrivals), seed=123, workers=2)
forest = assemble(plan, bank, workers=2)
```

The tree and Poisson streams are separate even with the same root seed. A larger
pre-generated iid bank is also valid; only its first N entries are consumed.
Do not pick trees by their survival outcomes, recycle entries, or assemble a bank
of fixed-root test trees under the uniform-immigration certificate. For paired
exit comparisons, make a new plan for the new model's certificate, extend the iid
bank if needed, and reuse trees only across, not within, populations. Such paired
outputs are correlated. Fragmentation parameters must stay unchanged in v0.1.

## Reproducibility and resource safety

Each family has streams indexed by (component, family ID, root seed), using
NumPy SeedSequence/PCG64. One and two workers, or different batch sizes, produce
identical outputs within the same numerical environment and execution mode.
Different seeds are needed for independent Monte Carlo replicates. Reproducibility
across NumPy versions or architectures is not promised.

Generation and pruning use batched ProcessPoolExecutor with explicit `spawn`.
Only a bounded number of batches is in flight. Defaults abort above 100000
immigrants, 1000000 nodes per family, or 2000000 stored nodes in a forest. These are
resource protections, not modelling thresholds; exceeding them raises an error
rather than returning a biased partial sample. Peak worker memory is additional
to the final stored-array size. Raising caps requires checking available memory.
Do not condition scientific results on which seeded runs fit a resource budget.

No rate clipping, artificial minimum exit, time discretisation, tree-recycling
approximation or background task is hidden in the implementation. Unrepresentable
clocks or split endpoints raise rather than silently altering the stochastic law.

## What is not yet implemented

GPU/compiled kernels, persistent on-disk tree banks, cross-size/time retiming,
clock-integrated stationary-mean estimation, arbitrary immigration distributions,
clinical fitting, and perfect infinite-past sampling. This is a tested reference
implementation before optimisation, not a production release or a claim of
measured Mac mini speedup.

Official implementation references:
- NumPy parallel random generation: https://numpy.org/doc/stable/reference/random/parallel.html
- Python process pools: https://docs.python.org/3/library/concurrent.futures.html

[Daily review handoff](docs/WORKFLOW.md) describes the small-task workflow. No
recurring task or automatic repository write has been activated by this package.
