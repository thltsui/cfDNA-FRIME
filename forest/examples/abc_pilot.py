"""Bounded rejection ABC for separate mtDNA length profiles; see docs/ABC_PILOT.md.

Stages: probe (engineering only), infer (synthetic + calibration + rejection),
ppc (fresh posterior-predictive populations). No failed realization is retried.
"""
from __future__ import annotations
import argparse
import csv
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
from time import perf_counter

import numpy as np
from frime_forest import Exit, Model, MonteCarlo, replicate_seed
from frime_forest.abc import observe_counts, cdf_distance, rejection_indices
from frime_forest.certificate import certify
from frime_forest.trees import rng_for
from monte_carlo import provenance

L_BP = 16569
CUTOFF_BP = 39.5
NAMES = ('alpha', 'boundary_bp', 'exit_at_half_boundary', 'immigration')
TRUTH = np.array([1., 500., 1.5, 30.])
PRIOR = dict(alpha=['uniform', 0, 2], boundary_bp=['log_uniform', 50, 6000],
             exit_at_half_boundary=['log_uniform', .5, 4], immigration=['log_uniform', 10, 100])
STAGES = dict(probe=100, truth=200, calibration=300, production=400, ppc=500)


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def edges_bp():
    # Recorded integer lengths: <=k corresponds to a continuous upper edge k+.5.
    upper = np.unique(np.r_[np.rint(np.geomspace(40, L_BP, 49)).astype(int),
                            [100, 200, 500, 1000, 2000, 5000, 10000]])
    return np.r_[CUTOFF_BP, upper + .5]


def model_for(theta, *, cutoff_bp=CUTOFF_BP):
    alpha, boundary_bp, e, nu = map(float, theta)
    b = boundary_bp / L_BP
    return Model(L=1, ell=cutoff_bp/L_BP, c_f=b**(-alpha), alpha=alpha,
                 a=1, b=1, c_i=nu, exit=Exit('pfb', c_e=e*b, power=-1, boundary=b))


def prior_draws(n, seed, stage):
    rng = rng_for(seed, stage, 10)
    return np.column_stack([rng.uniform(0, 2, n), np.exp(rng.uniform(np.log(50), np.log(6000), n)),
                            np.exp(rng.uniform(np.log(.5), np.log(4), n)),
                            np.exp(rng.uniform(np.log(10), np.log(100), n))])


def read_data(path, edges):
    groups = {}
    keys = set()
    with path.open(newline='') as f:
        for row in csv.DictReader(f):
            if row['chr'] != 'chrM':
                continue
            x, w = int(row['length']), int(row['freq'])
            if x <= 0 or w <= 0 or (row['id'], x) in keys:
                raise ValueError('invalid or duplicate mtDNA frequency row')
            keys.add((row['id'], x))
            groups.setdefault(row['id'], []).append((x, w))
    records, counts = [], []
    for i, sample in enumerate(sorted(groups)):
        xy = np.array(groups[sample], dtype=np.int64)
        h = np.histogram(xy[:, 0], edges, weights=xy[:, 1])[0]
        if xy[:, 0].max() > L_BP:
            raise ValueError('observed length exceeds declared immigrant maximum')
        if not h.sum():
            raise ValueError('sample has no reads in the analysis window')
        records.append(dict(label=f'S{i+1:02d}', source_id=sample, n=int(h.sum()),
                            source_n=int(xy[:, 1].sum()), excluded=int(xy[:, 1].sum()-h.sum())))
        counts.append(h)
    if not counts:
        raise ValueError('no mtDNA profiles found')
    return records, np.array(counts)


def source_info():
    info = provenance()
    info['pilot_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return info


def preflight(theta, delta, cap):
    certificates = [certify(model_for(t), delta) for t in theta]
    expected = np.array([c.model.c_i*c.T for c in certificates])
    if np.any(expected > cap):
        raise RuntimeError('prior proposal exceeds family cap; whole experiment aborted, no filtering')
    return dict(proposals=len(theta), expected_families_quantiles=np.quantile(expected, [0, .5, .95, 1]).tolist(),
                T_quantiles=np.quantile([c.T for c in certificates], [0, .5, .95, 1]).tolist())


class RunLog:
    def __init__(self, runner, args, log):
        self.runner, self.args, self.log = runner, args, log
        self.n = 0
        self.nodes = 0
        self.immigrants = 0
        self.lambda_sum = 0.

    def population(self, theta, stage, index):
        root = replicate_seed(self.args.seed, stage)
        result = self.runner.run(model_for(theta), 1, bins=edges_bp()/L_BP,
                                 seed=root, start_replicate=index, delta=self.args.delta)
        h = result.histograms[0]
        self.n += 1
        self.nodes += int(result.generated_nodes[0])
        self.immigrants += int(result.immigrants[0])
        self.lambda_sum += result.certificate.lambda_bound
        self.log.write(json.dumps(dict(manifest=('ppc-manifest.json' if self.args.stage == 'ppc' else 'manifest.json'), stage=stage, proposal=index,
            theta=dict(zip(NAMES, map(float, theta))), started_utc=result.started_utc,
            seed=str(result.seeds[0]), ensemble_seed=str(root), workers=self.args.workers,
            batch_size=self.args.batch_size, mode='lazy', model=asdict(result.model),
            certificate=asdict(result.certificate), histogram=h.tolist(),
            observable_count=int(result.counts[0]), observable_mass=float(result.masses[0]),
            immigrants=int(result.immigrants[0]), generated_nodes=int(result.generated_nodes[0]),
            seconds=float(result.seconds[0])), allow_nan=False) + '\n')
        self.log.flush()
        return h

    def totals(self):
        return dict(populations=self.n, immigrants=self.immigrants, generated_nodes=self.nodes,
                    history_union_bound=min(1., self.lambda_sum))


def bank(run, theta, stage, observed):
    hist = np.empty((len(theta), observed.shape[1]), dtype=np.int64)
    distances = np.empty((len(theta), len(observed)))
    for i, t in enumerate(theta):
        h = hist[i] = run.population(t, stage, i)
        for j, y in enumerate(observed):
            z = observe_counts(h, int(y.sum()), rng_for(run.args.seed, i, stage*1000+j))
            distances[i, j] = cdf_distance(z, y)
        if (i+1) % 32 == 0:
            print(json.dumps(dict(stage=stage, completed=i+1, total=len(theta))), flush=True)
    return hist, distances


def wilson(k, n):
    z = 1.959963984540054
    p = k/n
    center = (p+z*z/(2*n))/(1+z*z/n)
    radius = z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/(1+z*z/n)
    return [float(center-radius), float(center+radius)]


def infer(run, args, output):
    labels, observed = read_data(args.data, edges_bp())
    truth_h = run.population(TRUTH, STAGES['truth'], 0)
    truth_y = observe_counts(truth_h, 1000, rng_for(args.seed, 0, 200000))
    if not truth_y.sum():
        raise RuntimeError('synthetic truth produced no reads; abort diagnostic, do not retry')
    labels.append(dict(label='synthetic', n=1000, truth=dict(zip(NAMES, TRUTH.tolist()))))
    observed = np.vstack([observed, truth_y])
    dump(output/'datasets.json', labels)
    np.savez_compressed(output/'observed.npz', bins_bp=edges_bp(), counts=observed,
                        labels=np.array([x['label'] for x in labels]))
    calibration = prior_draws(args.calibration, args.seed, STAGES['calibration'])
    production = prior_draws(args.proposals, args.seed, STAGES['production'])
    dump(output/'preflight.json', {s:preflight(t, args.delta, args.max_families)
         for s,t in [('calibration',calibration), ('production',production)]})
    h_cal, d_cal = bank(run, calibration, STAGES['calibration'], observed)
    # Freeze an order statistic before simulating any production proposal.
    rank = math.ceil(.05*args.calibration)-1
    epsilon = np.sort(d_cal, axis=0)[rank]
    if not np.all(np.isfinite(epsilon)):
        raise RuntimeError('calibration tolerance is not finite; no production inference')
    dump(output/'tolerances.json', dict(rule='ceil(0.05*n)th calibration distance; frozen before production',
                                       epsilon=epsilon.tolist()))
    np.savez_compressed(output/'calibration.npz', theta=calibration, population_histograms=h_cal, distances=d_cal)
    h_prod, d_prod = bank(run, production, STAGES['production'], observed)
    np.savez_compressed(output/'production.npz', theta=production, population_histograms=h_prod,
                        distances=d_prod, epsilon=epsilon)
    summaries = []
    for j, label in enumerate(labels):
        ix = rejection_indices(d_prod[:,j], epsilon[j])
        quantiles = np.quantile(production[ix], [.025, .5, .975], axis=0).tolist() if len(ix) else None
        summary = dict(label=label['label'], n=label['n'], epsilon=float(epsilon[j]),
                       accepted=len(ix), proposed=args.proposals,
                       acceptance_wilson95=wilson(len(ix), args.proposals),
                       minimum_distance=(float(np.min(d_prod[:,j])) if np.any(np.isfinite(d_prod[:,j])) else None),
                       no_read_outcomes=int(np.sum(~np.isfinite(d_prod[:,j]))),
                       distance_le_005=int(np.sum(d_prod[:,j] <= .05)),
                       distance_le_010=int(np.sum(d_prod[:,j] <= .10)),
                       parameter_names=NAMES, parameter_quantiles_025_50_975=quantiles,
                       accepted_proposal_ids=ix.tolist())
        if label['label'] == 'synthetic' and len(ix):
            q=np.array(quantiles)
            summary['truth_in_marginal_95_intervals'] = ((q[0] <= TRUTH)&(TRUTH <= q[2])).tolist()
            transformed=q.copy(); transformed[:,1:]=np.log(transformed[:,1:])
            prior95_width=.95*np.array([2., np.log(6000/50), np.log(4/.5), np.log(100/10)])
            summary['posterior_to_prior_95_width_ratio_uniform_or_log_scale'] = ((transformed[2]-transformed[0])/prior95_width).tolist()
        summaries.append(summary)
    dump(output/'posterior.json', summaries)


def ppc(run, args, output):
    data = np.load(output/'observed.npz')
    prod = np.load(output/'production.npz')
    summaries = json.loads((output/'posterior.json').read_text())
    observations = data['counts']
    results = []
    for j, summary in enumerate(summaries):
        ix = np.array(summary['accepted_proposal_ids'], dtype=int)
        if len(ix) == 0:
            results.append(dict(label=summary['label'], status='no posterior particles'))
            continue
        chosen = rng_for(args.seed, j, 501000).choice(ix, size=args.ppc, replace=True)
        ys, hs = [], []
        stage = STAGES['ppc']+j
        for i, source_index in enumerate(chosen):
            h = run.population(prod['theta'][source_index], stage, i)
            y = observe_counts(h, int(observations[j].sum()), rng_for(args.seed, i, stage*1000+999))
            hs.append(h); ys.append(y)
        ys, hs = np.array(ys), np.array(hs)
        good = ys.sum(axis=1)>0
        # No-read outcomes stay in stored output and unconditional failure counts.
        cdfs = np.full(ys.shape, np.nan)
        cdfs[good] = np.cumsum(ys[good],axis=1)/ys[good].sum(axis=1)[:,None]
        if not np.any(good):
            np.savez_compressed(output/f'ppc-{summary["label"]}.npz', counts=ys,
                population_histograms=hs, production_proposal_ids=chosen)
            results.append(dict(label=summary['label'], status='all posterior predictive outcomes have no reads',
                                predictive_draws=args.ppc, no_read_outcomes=args.ppc))
            continue
        center = np.mean(cdfs[good], axis=0)
        radial = np.max(np.abs(cdfs[good]-center), axis=1)
        radius = float(np.quantile(radial, .95, method='higher'))
        obs_cdf = np.cumsum(observations[j])/observations[j].sum()
        obs_gap = float(np.max(np.abs(obs_cdf-center)))
        rank = (1+np.sum(radial >= obs_gap))/(1+len(radial))
        direct = np.array([cdf_distance(y,observations[j]) for y in ys])
        finite = direct[np.isfinite(direct)]
        result = dict(label=summary['label'], predictive_draws=args.ppc, no_read_outcomes=int((~good).sum()),
                      cdf_center_max_gap=obs_gap, simultaneous95_radius=radius,
                      observed_inside_simultaneous95=bool(obs_gap<=radius),
                      descriptive_predictive_tail_rank=float(rank),
                      direct_distance_quantiles_05_50_95=np.quantile(finite,[.05,.5,.95]).tolist(),
                      note='Descriptive in-sample posterior predictive check; band/rank conditional on nonempty outcome, not a calibrated p-value or external validation.')
        results.append(result)
        np.savez_compressed(output/f'ppc-{summary["label"]}.npz', counts=ys, population_histograms=hs,
            production_proposal_ids=chosen, center_cdf=center, radius=radius,
            pointwise025=np.quantile(cdfs[good],.025,axis=0),pointwise975=np.quantile(cdfs[good],.975,axis=0))
        print(json.dumps(result), flush=True)
    dump(output/'predictive.json', results)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage',choices=['probe','infer','ppc'],required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--data',type=Path,default=Path(__file__).resolve().parents[2]/'data_preprocessing/mtDNA_modeling_data.csv')
    p.add_argument('--seed',type=int,default=20261012)
    p.add_argument('--workers',type=int,default=4)
    p.add_argument('--batch-size',type=int,default=128)
    p.add_argument('--delta',type=float,default=1e-10)
    p.add_argument('--probe-size',type=int,default=32)
    p.add_argument('--calibration',type=int,default=256)
    p.add_argument('--proposals',type=int,default=2048)
    p.add_argument('--ppc',type=int,default=100)
    p.add_argument('--max-families',type=int,default=100000)
    p.add_argument('--max-total-nodes',type=int,default=10000000)
    args=p.parse_args()
    for name in ('probe_size','calibration','proposals','ppc','workers','batch_size'):
        if getattr(args,name)<1:
            p.error(f'{name} must be positive')
    output=args.output
    if args.stage=='ppc':
        manifest=json.loads((output/'manifest.json').read_text())
        if not (output/'infer-complete.json').exists() or (output/'infer-FAILED.json').exists():
            p.error('PPC requires a completed, successful inference stage')
        current = source_info()
        for key in ('source_sha256', 'pilot_sha256'):
            if current[key] != manifest['provenance'][key]:
                p.error('Scientific source changed since inference; start a new experiment')
        if not np.array_equal(np.load(output/'observed.npz')['bins_bp'], edges_bp()):
            p.error('Observation bins changed since inference')
        old=manifest['arguments']
        for name in ('seed','delta','max_families','max_total_nodes'):
            if getattr(args,name)!=old[name]:
                p.error(f'PPC {name} differs from inference')
        if (output/'predictive.json').exists() or list(output.glob('ppc-S*.npz')):
            p.error('PPC output already exists; refusing overwrite/retry of partial run')
        dump(output/'ppc-manifest.json',dict(provenance=source_info(),arguments=vars(args)|dict(output=str(output),data=str(args.data))))
    else:
        output.mkdir(parents=True,exist_ok=False)
        manifest=dict(provenance=source_info(),timestamp_utc=datetime.now(timezone.utc).isoformat(),
            data_sha256=hashlib.sha256(args.data.read_bytes()).hexdigest(), prior=PRIOR,
            fixed=dict(L_bp=L_BP,biological_cutoff_bp=1,projected_cutoff_bp=CUTOFF_BP,
                       split_a=1,split_b=1,exit='PFB',exit_power=-1),
            observation='Multinomial count weighting within rounded lengths >=40bp; empty population yields no-read outcome',
            discrepancy='maximum binned CDF gap', calibration_quantile=.05,
            notes=['Separate marginal ABC fits share prior-predictive bank; cross-profile Monte Carlo is correlated.',
                   'rho(B)=1 fixes time units. Immigration is total intensity, not per-unit-length.',
                   'Exact downward-fragmentation projection; certificate covers observable population only.',
                   'Calibration quantile controls acceptance, not goodness of fit. Synthetic check is not SBC.',
                   'Resource/numerical failures abort the stage, never become ABC rejections.'],
            arguments=vars(args)|dict(output=str(output),data=str(args.data)))
        dump(output/'manifest.json',manifest)
    started=perf_counter()
    try:
        with (output/f'{args.stage}-runs.jsonl').open('x') as log, MonteCarlo(workers=args.workers,
            batch_size=args.batch_size,max_families=args.max_families,max_total_nodes=args.max_total_nodes) as runner:
            run=RunLog(runner,args,log)
            if args.stage=='probe':
                theta=prior_draws(args.probe_size,args.seed,STAGES['probe'])
                dump(output/'preflight.json',preflight(theta,args.delta,args.max_families))
                for i,t in enumerate(theta):
                    run.population(t,STAGES['probe'],i)
                    if (i+1)%8==0:
                        print(f'engineering probe {i+1}/{len(theta)}',flush=True)
            elif args.stage=='infer':
                infer(run,args,output)
            else:
                ppc(run,args,output)
            totals=run.totals()
        totals['wall_seconds']=perf_counter()-started
        if args.stage == 'infer':
            datasets=len(json.loads((output/'datasets.json').read_text()))
            totals['planned_whole_pipeline_history_union_bound'] = min(1.,args.delta*(1+args.calibration+args.proposals+datasets*args.ppc))
            totals['planned_ppc_per_dataset']=args.ppc
        elif args.stage == 'ppc':
            datasets=len(json.loads((output/'datasets.json').read_text()))
            totals['whole_pipeline_history_union_bound'] = min(1.,args.delta*(1+old['calibration']+old['proposals']+datasets*args.ppc))
        dump(output/f'{args.stage}-complete.json',totals)
        print(json.dumps(totals),flush=True)
    except BaseException as exc:
        dump(output/f'{args.stage}-FAILED.json',dict(error=repr(exc),elapsed_seconds=perf_counter()-started,
             note='Partial output is diagnostic only; no completed scientific result. No seed was skipped.'))
        raise

if __name__=='__main__':
    main()
