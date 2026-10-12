"""Execute a frozen four-case immigration sensitivity plan; not model evidence.

Each plan row specifies one full Model and its parameter coordinates. Screening
uses one independent random population/read outcome per row. Confirmation uses
100 fresh populations at each of three best screened settings per case.
"""
from __future__ import annotations

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from frime_forest import Exit, Immigration, Model, MonteCarlo, replicate_seed
from frime_forest.abc import cdf_distance, observe_counts
from frime_forest.certificate import certify
from frime_forest.trees import rng_for
from abc_pilot import CUTOFF_BP, L_BP, dump, edges_bp, read_data, wilson
from monte_carlo import provenance

THRESHOLDS = (.05, .10, .15)
TAIL_BP = (1000, 5000, 10000)


def model_from(fields):
    fields = dict(fields)
    fields['exit'] = Exit(**fields['exit'])
    fields['immigration'] = Immigration(**fields.get('immigration', {}))
    return Model(**fields)


def source_info():
    info = provenance()
    for name in ('immigration_sensitivity.py', 'plan_immigration_sensitivity.py', 'abc_pilot.py', 'monte_carlo.py'):
        p = Path(__file__).with_name(name)
        info['source_sha256']['forest/examples/'+name] = hashlib.sha256(p.read_bytes()).hexdigest()
    return info


def hashes(output, names):
    return {name:hashlib.sha256((output/name).read_bytes()).hexdigest() for name in names}


def select_candidates(rows, distances, top):
    selected=[]
    for case in dict.fromkeys(row['case'] for row in rows):
        ix=[i for i,row in enumerate(rows) if row['case']==case]
        order=sorted(ix,key=lambda i:(distances[i],rows[i]['id']))
        selected.extend(rows[i] for i in [i for i in order if np.isfinite(distances[i])][:top])
    return selected


def success_summary(d):
    d = np.asarray(d)
    finite = d[np.isfinite(d)]
    return dict(n=len(d), no_read_outcomes=int(np.sum(~np.isfinite(d))),
                distance_quantiles_05_50_95=(np.quantile(finite, [.05,.5,.95]).tolist() if len(finite) else None),
                minimum_distance=(float(finite.min()) if len(finite) else None),
                successes={str(t):dict(count=int(np.sum(d<=t)),
                    proportion=float(np.mean(d<=t)), wilson95=wilson(int(np.sum(d<=t)),len(d)))
                    for t in THRESHOLDS})


def predictive_summary(y, observed):
    good = y.sum(axis=1)>0
    d = np.array([cdf_distance(row,observed) for row in y])
    result = success_summary(d)
    if np.any(good):
        cdf = np.cumsum(y[good],axis=1)/y[good].sum(axis=1)[:,None]
        observed_cdf = observed.cumsum()/observed.sum()
        mean = cdf.mean(axis=0)
        se = cdf.std(axis=0,ddof=1)/np.sqrt(len(cdf)) if len(cdf)>1 else np.full(len(mean),np.nan)
        gap = np.abs(mean-observed_cdf)
        ix = int(np.argmax(gap))
        result.update(mean_cdf=mean.tolist(), mean_cdf_max_gap=float(gap[ix]),
                      gap_at_bp=float(edges_bp()[1:][ix]-.5),
                      cdf_mean_mcse_at_max_gap=(float(se[ix]) if len(cdf)>1 else None),
                      pointwise025=np.quantile(cdf,.025,axis=0).tolist(),
                      pointwise975=np.quantile(cdf,.975,axis=0).tolist())
        tails = {}
        for threshold in TAIL_BP:
            j = int(np.flatnonzero(edges_bp()[1:]==threshold+.5)[0])
            values = 1-cdf[:,j]
            tails[str(threshold)] = dict(observed=float(1-observed_cdf[j]),
                predictive_mean=float(values.mean()), mcse=(float(se[j]) if len(cdf)>1 else None),
                predictive025_975=np.quantile(values,[.025,.975]).tolist())
        result['tails'] = tails
    result['interpretation'] = ('Fresh in-sample predictive diagnostic at selected fixed parameters; '
        'not a calibrated p-value, posterior model probability, or external validation. '
        'CDF summaries condition on nonempty reads; unconditional no-read counts retained.')
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage',choices=('check','screen','confirm'),required=True)
    p.add_argument('--plan',type=Path,help='Required only for check; copied and frozen in output')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--data',type=Path,default=Path(__file__).resolve().parents[2]/'data_preprocessing/mtDNA_modeling_data.csv')
    p.add_argument('--workers',type=int,default=8)
    p.add_argument('--batch-size',type=int,default=128)
    args = p.parse_args()
    output = args.output
    if args.stage=='check':
        if args.plan is None:
            p.error('check requires --plan')
        output.mkdir(parents=True,exist_ok=False)
        plan = json.loads(args.plan.read_text())
        dump(output/'plan.json',plan)
        labels, counts = read_data(args.data,edges_bp())
        target = next(i for i,r in enumerate(labels) if r['label']==plan['target'])
        observed = counts[target]
        dump(output/'dataset.json',labels[target])
        np.savez_compressed(output/'observed.npz',counts=observed,bins_bp=edges_bp())
        manifest = dict(provenance=source_info(),timestamp_utc=datetime.now(timezone.utc).isoformat(),
            data_sha256=hashlib.sha256(args.data.read_bytes()).hexdigest(),
            plan_sha256=hashlib.sha256((output/'plan.json').read_bytes()).hexdigest(),
            observation_sha256=hashes(output,('observed.npz','dataset.json')),
            workers=args.workers,batch_size=args.batch_size,
            note='Frozen exploratory comparison, not a posterior model-selection procedure.')
        dump(output/'manifest.json',manifest)
    else:
        plan = json.loads((output/'plan.json').read_text())
        manifest = json.loads((output/'manifest.json').read_text())
        if (source_info()['source_sha256'] != manifest['provenance']['source_sha256'] or
            hashlib.sha256(args.data.read_bytes()).hexdigest()!=manifest['data_sha256'] or
            hashlib.sha256((output/'plan.json').read_bytes()).hexdigest()!=manifest['plan_sha256']):
            p.error('scientific source, data or plan changed; refuse mixed experiment')
        if hashes(output,manifest['observation_sha256'])!=manifest['observation_sha256']:
            p.error('frozen observation changed')
        if not (output/'check-complete.json').exists():
            p.error('successful preflight required')
        observed = np.load(output/'observed.npz')['counts']
    rows = plan['rows']
    for key in ('confirmation_repeats','top_per_case','max_families','max_nodes','max_total_nodes'):
        if type(plan[key]) is not int or plan[key]<1:
            p.error(f'{key} must be a positive integer')
    if not rows:
        p.error('nonempty plan required')
    cases = list(dict.fromkeys(row['case'] for row in rows))
    if len({row['id'] for row in rows})!=len(rows):
        p.error('plan IDs must be unique')
    if any(model_from(row['model']).L != 1 or
           model_from(row['model']).ell != CUTOFF_BP/L_BP for row in rows):
        p.error('plan must use the declared observable projection')
    for row in rows:
        if type(row['id']) is not int or row['id']<0:
            p.error('nonnegative integer plan IDs required')
        if (row['id']+1)*plan['confirmation_repeats']>2**63:
            p.error('planned replicate IDs must fit signed int64')
    if args.stage=='confirm' and not (output/'screen-complete.json').exists():
        p.error('confirmation requires successful complete screening')
    if args.stage=='confirm':
        done=json.loads((output/'screen-complete.json').read_text())
        if hashes(output,done['artifacts_sha256'])!=done['artifacts_sha256']:
            p.error('screening artifacts changed')
        saved=json.loads((output/'selected.json').read_text())
        distances=np.load(output/'screen.npz')['distances']
        if saved!=select_candidates(rows,distances,plan['top_per_case']):
            p.error('selected settings do not match the frozen screening rule')
    if (output/f'{args.stage}-started.json').exists():
        p.error('stage already attempted; preserve outputs, never silently retry')
    dump(output/f'{args.stage}-started.json',dict(provenance=source_info(),workers=args.workers,
         batch_size=args.batch_size,timestamp_utc=datetime.now(timezone.utc).isoformat()))
    start = perf_counter()
    totals = dict(populations=0,immigrants=0,generated_nodes=0,history_union_bound=0.)
    try:
        if args.stage=='check':
            records=[]
            for row in rows:
                m=model_from(row['model']); c=certify(m,plan['delta'])
                records.append(dict(id=row['id'],case=row['case'],expected_families=m.c_i*c.T,
                                    certificate=asdict(c)))
            dump(output/'certificates.json',records)
            summary={case:dict(proposals=sum(r['case']==case for r in records),
                expected_families_quantiles=np.quantile([r['expected_families'] for r in records if r['case']==case],[0,.5,.95,1]).tolist(),
                total_expected_families=sum(r['expected_families'] for r in records if r['case']==case),
                exceeding_cap=sum(r['expected_families']>plan['max_families'] for r in records if r['case']==case)) for case in cases}
            dump(output/'preflight.json',summary)
            print(json.dumps(summary),flush=True)
            if any(v['exceeding_cap'] for v in summary.values()):
                raise RuntimeError('whole plan fails prospective resource check; no proposals filtered or simulated')
        else:
            selected = (rows if args.stage=='screen' else json.loads((output/'selected.json').read_text()))
            nrepeat = 1 if args.stage=='screen' else plan['confirmation_repeats']
            stage = 710 if args.stage=='screen' else 720
            ensemble_seed = replicate_seed(plan['seed'],stage)
            hs=[]; ys=[]; distances=[]; output_rows=[]; summaries=[]
            with (output/f'{args.stage}-runs.jsonl').open('x') as log, MonteCarlo(
                workers=args.workers,batch_size=args.batch_size,max_families=plan['max_families'],
                max_nodes=plan['max_nodes'],max_total_nodes=plan['max_total_nodes']) as mc:
                for k,row in enumerate(selected):
                    local_y=[]
                    for j in range(nrepeat):
                        index=row['id']*plan['confirmation_repeats']+j
                        m=model_from(row['model'])
                        result=mc.run(m,1,bins=edges_bp()/L_BP,seed=ensemble_seed,
                                      start_replicate=index,delta=plan['delta'])
                        h=result.histograms[0]
                        y=observe_counts(h,int(observed.sum()),rng_for(plan['seed'],index,stage))
                        distance=cdf_distance(y,observed)
                        record=dict(case=row['case'],proposal=row['id'],repeat=j,theta=row['theta'],
                            model=asdict(m),certificate=asdict(result.certificate),
                            started_utc=result.started_utc,seed=str(result.seeds[0]),
                            ensemble_seed=str(ensemble_seed),observation_stream=[plan['seed'],index,stage],
                            workers=args.workers,batch_size=args.batch_size,mode='lazy',
                            histogram=h.tolist(),reads=y.tolist(),
                            distance=(distance if np.isfinite(distance) else None),
                            immigrants=int(result.immigrants[0]),generated_nodes=int(result.generated_nodes[0]),
                            observable_count=int(result.counts[0]),observable_mass=float(result.masses[0]),
                            seconds=float(result.seconds[0]))
                        log.write(json.dumps(record,allow_nan=False)+'\n');log.flush()
                        hs.append(h);ys.append(y);local_y.append(y);distances.append(distance)
                        output_rows.append((row['id'],j))
                        totals['populations']+=1
                        totals['immigrants']+=int(result.immigrants[0])
                        totals['generated_nodes']+=int(result.generated_nodes[0])
                        totals['history_union_bound']+=result.certificate.lambda_bound
                    if args.stage=='confirm':
                        s=predictive_summary(np.array(local_y),observed)
                        s.update(case=row['case'],proposal=row['id'],theta=row['theta'])
                        summaries.append(s)
                        dump(output/'confirmation-partial.json',summaries)
                    if (k+1)%16==0 or args.stage=='confirm':
                        print(json.dumps(dict(stage=args.stage,completed=k+1,total=len(selected),
                                             wall_seconds=perf_counter()-start)),flush=True)
            np.savez_compressed(output/f'{args.stage}.npz',histograms=np.array(hs),reads=np.array(ys),
                                distances=np.array(distances),indices=np.array(output_rows))
            if args.stage=='screen':
                summaries=[]
                for case in cases:
                    ix=[i for i,row in enumerate(rows) if row['case']==case]
                    s=success_summary(np.array(distances)[ix]);s['case']=case;summaries.append(s)
                selected=select_candidates(rows,distances,plan['top_per_case'])
                dump(output/'selected.json',selected)
                dump(output/'screen-summary.json',summaries)
                totals['artifacts_sha256']=hashes(output,('selected.json','screen.npz','screen-summary.json'))
            else:
                dump(output/'confirmation-summary.json',summaries)
        totals['wall_seconds']=perf_counter()-start
        dump(output/f'{args.stage}-complete.json',totals)
        print(json.dumps(totals),flush=True)
    except BaseException as exc:
        dump(output/f'{args.stage}-FAILED.json',dict(error=repr(exc),**totals,
             wall_seconds=perf_counter()-start,note='Partial output diagnostic only; no seeds skipped.'))
        raise


if __name__=='__main__':
    main()
