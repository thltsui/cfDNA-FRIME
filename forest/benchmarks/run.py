"""Reproducible cold-start benchmark, including spawn and IPC (not GPU timings).

From forest/: python benchmarks/run.py --workers 1 2 4 --repeats 3
The report describes its actual host; no Mac mini extrapolation is performed.
"""
from __future__ import annotations
import argparse
from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
import platform
import statistics
import sys
from time import perf_counter

import numpy as np
from frime_forest import Exit, Model, sample


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", nargs="+", type=int, default=[1,2])
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--immigration", type=float, default=100)
    parser.add_argument("--cutoff", type=float, default=.02)
    parser.add_argument("--output", type=Path, default=Path("results/benchmark.json"))
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("repeats must be positive")
    model = Model(c_i=args.immigration,ell=args.cutoff,alpha=1,exit=Exit("pfb",power=-1,boundary=.4))
    records = []
    fingerprints = {}
    for mode in ("lazy","skeleton"):
        for workers in args.workers:
            timings, runs = [], []
            for repeat in range(args.repeats):
                seed = 12345+repeat
                start = perf_counter()
                forest = sample(model,seed=seed,workers=workers,mode=mode,delta=1e-6)
                elapsed = perf_counter()-start
                values = forest.at()
                digest = hashlib.sha256(values.tobytes()).hexdigest()
                key = mode,seed
                if key in fingerprints and digest != fingerprints[key]:
                    raise AssertionError("worker count changed the sample")
                fingerprints[key] = digest
                timings.append(elapsed)
                runs.append(dict(seed=seed,seconds=elapsed,immigrants=len(forest.families),
                                 stored_nodes=forest.stored_nodes,fragments=len(values),sha256=digest))
            record = dict(mode=mode,workers=workers,median_seconds=statistics.median(timings),runs=runs)
            records.append(record)
            print(f"{mode:8s} workers={workers}: median {record['median_seconds']:.3f} s",flush=True)
    report = dict(model=asdict(model),certificate=asdict(forest.plan.certificate),
                  host=dict(platform=platform.platform(),machine=platform.machine(),cpu_count=os.cpu_count(),
                            python=sys.version.split()[0],numpy=np.__version__),
                  includes="cold process startup, IPC, generation, pruning, assembly; excludes snapshot query",
                  legacy_comparison=False,records=records)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+"\n")


if __name__ == "__main__":
    main()
