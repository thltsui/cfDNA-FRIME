"""Run from forest/: python examples/demo.py --workers 2 --mode lazy."""
from __future__ import annotations
import argparse
from dataclasses import asdict
import json
import platform
from pathlib import Path
import sys
from time import perf_counter

import numpy as np
from frime_forest import Exit, Model, sample


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--mode", choices=("lazy", "skeleton"), default="lazy")
    parser.add_argument("--exit", choices=("constant", "pnb", "pfb"), default="constant")
    parser.add_argument("--alpha", type=float, default=0.0)
    parser.add_argument("--cutoff", type=float, default=0.01)
    parser.add_argument("--immigration", type=float, default=20.0)
    parser.add_argument("--delta", type=float, default=1e-6)
    parser.add_argument("--horizon", type=float, default=2.0)
    parser.add_argument("--seed", type=int, default=20260908)
    parser.add_argument("--output", type=Path, default=Path("results/demo"))
    args = parser.parse_args()
    model = Model(ell=args.cutoff, alpha=args.alpha, c_i=args.immigration, exit=Exit(args.exit))
    start = perf_counter()
    forest = sample(model, delta=args.delta, horizon=args.horizon, seed=args.seed,
                    workers=args.workers, mode=args.mode)
    elapsed = perf_counter()-start
    times = np.linspace(0, args.horizon, 5)
    bins = np.geomspace(model.ell, model.L, 41)
    histograms = np.array([forest.histogram(bins, float(t)) for t in times])
    snapshot = forest.at(0)
    metadata = dict(model=asdict(model), certificate=asdict(forest.plan.certificate),
                    seed=args.seed, workers=args.workers, mode=args.mode,
                    elapsed_seconds=elapsed, immigrants=len(forest.families),
                    stored_nodes=forest.stored_nodes, snapshot_count=len(snapshot),
                    snapshot_mass=float(snapshot.sum()), platform=platform.platform(),
                    python=sys.version.split()[0], numpy=np.__version__)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/"metadata.json").write_text(json.dumps(metadata, indent=2)+"\n")
    np.savez_compressed(args.output/"snapshot.npz", sizes=snapshot, times=times,
                        bins=bins, histograms=histograms, arrivals=forest.plan.arrivals)
    print(json.dumps(metadata, indent=2))


if __name__ == "__main__":
    main()
