"""Save independent FRIME Monte Carlo populations, summaries and provenance.

From forest/: python examples/monte_carlo.py --n 1000 --workers 4 --output results/mc-001
Defaults are dimensionless illustrative parameters, not an mtDNA calibration.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
from time import perf_counter

import numpy as np
from frime_forest import Exit, Immigration, Model, MonteCarlo


def command(args, cwd=None):
    try:
        return subprocess.check_output(args, cwd=cwd, stderr=subprocess.DEVNULL, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def provenance():
    root = Path(__file__).resolve().parents[2]
    return dict(git_commit=command(['git','rev-parse','HEAD'], root),
                git_status=command(['git','status','--porcelain'], root),
                source_sha256={str(f.relative_to(root)):hashlib.sha256(f.read_bytes()).hexdigest()
                               for f in sorted((root/'forest/frime_forest').glob('*.py'))},
                machine=dict(platform=platform.platform(), architecture=platform.machine(),
                             cpu_count=os.cpu_count(),
                             chip=command(['sysctl','-n','machdep.cpu.brand_string']),
                             ram_bytes=command(['sysctl','-n','hw.memsize'])),
                python=platform.python_version(), numpy=np.__version__)


def save_result(result, output, metadata):
    # A fresh directory prevents accidental replacement of a previous experiment.
    output.mkdir(parents=True, exist_ok=False)
    np.savez_compressed(output/'populations.npz', bins=result.bins,
        replicate_ids=result.replicate_ids, seeds=np.array([str(s) for s in result.seeds]),
        histograms=result.histograms, counts=result.counts, masses=result.masses,
        immigrants=result.immigrants, generated_nodes=result.generated_nodes,
        seconds=result.seconds)
    record = dict(**metadata, model=asdict(result.model), certificate=asdict(result.certificate),
                  started_utc=result.started_utc, seed=result.seed,
                  start_replicate=int(result.replicate_ids[0]), n=len(result.seeds),
                  workers=result.workers, batch_size=result.batch_size, mode=result.mode,
                  elapsed_seconds=result.elapsed_seconds,
                  per_replicate_delta=result.certificate.delta,
                  joint_history_bound=result.joint_history_bound,
                  statistics=result.statistics(), output='populations.npz',
                  notes=['Raw population counts, not normalized sequencing probabilities.',
                         'Caps are per population; failed ensembles are never returned as partial samples.',
                         'Monte Carlo SE is across independent populations, separate from history error.'])
    (output/'metadata.json').write_text(json.dumps(record, indent=2, allow_nan=False)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--n', type=int, default=1000)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--batch-size', type=int, default=128)
    parser.add_argument('--seed', type=int, default=20260914)
    parser.add_argument('--start-replicate', type=int, default=0)
    parser.add_argument('--delta', type=float, default=1e-6, help='history tolerance per replicate')
    parser.add_argument('--mode', choices=['lazy','skeleton'], default='lazy')
    parser.add_argument('--model-json', type=Path, help='Model fields; exit is a nested object')
    parser.add_argument('--bins-json', type=Path, help='JSON array of bin edges in model length units')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error('output already exists; choose a fresh experiment directory')
    if args.model_json:
        fields = json.loads(args.model_json.read_text())
        fields['exit'] = Exit(**fields.get('exit', {}))
        fields['immigration'] = Immigration(**fields.get('immigration', {}))
        model = Model(**fields)
    else:
        model = Model(c_i=100, ell=.02, alpha=1, exit=Exit('pfb', power=-1, boundary=.4))
    bins = (json.loads(args.bins_json.read_text()) if args.bins_json
            else np.linspace(model.ell, model.L, 51))
    metadata = provenance()
    start = perf_counter()
    with MonteCarlo(workers=args.workers, batch_size=args.batch_size) as runner:
        result = runner.run(model, args.n, bins=bins, seed=args.seed,
                            start_replicate=args.start_replicate, delta=args.delta, mode=args.mode)
    metadata['wall_with_pool_shutdown_seconds'] = perf_counter()-start
    save_result(result, args.output, metadata)
    print(json.dumps(dict(populations=args.n, seconds=result.elapsed_seconds,
                          joint_history_bound=result.joint_history_bound,
                          count=result.statistics()['count'], output=str(args.output)), indent=2))


if __name__ == '__main__':
    main()
