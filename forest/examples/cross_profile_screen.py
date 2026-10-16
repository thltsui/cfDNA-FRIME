"""Marginal ABC screens using a complete, audited random-population bank.

The bank is shared across profiles; this is not a joint/hierarchical likelihood.
No posterior-selected subset of populations may be used as the input bank.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from frime_forest.abc import cdf_distance, observe_counts
from frime_forest.trees import rng_for
from abc_pilot import dump, edges_bp, read_data
from immigration_sensitivity import source_info, success_summary, select_candidates

THRESHOLDS = (.05, .10, .15)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_bank(directory, data):
    """Reject incomplete, changed, filtered or mismatched scientific input."""
    if (directory/'screen-FAILED.json').exists():
        raise ValueError('failed bank cannot be screened')
    plan=json.loads((directory/'plan.json').read_text())
    manifest=json.loads((directory/'manifest.json').read_text())
    done=json.loads((directory/'screen-complete.json').read_text())
    if sha(directory/'plan.json')!=manifest['plan_sha256'] or sha(data)!=manifest['data_sha256']:
        raise ValueError('bank plan or data changed')
    for name,expected in {**manifest['observation_sha256'],**done['artifacts_sha256']}.items():
        if sha(directory/name)!=expected:
            raise ValueError('bank artifact changed: '+name)
    if source_info()['source_sha256']!=manifest['provenance']['source_sha256']:
        raise ValueError('scientific bank sources changed')
    with np.load(directory/'screen.npz',allow_pickle=False) as z:
        arrays={k:z[k].copy() for k in ('histograms','reads','distances','indices')}
    rows=plan['rows']; n=len(rows)
    if done['populations']!=n or arrays['indices'].tolist()!=[[r['id'],0] for r in rows]:
        raise ValueError('bank is not the complete original prior plan')
    if len({r['id'] for r in rows})!=n:
        raise ValueError('duplicate proposal IDs')
    observations=np.load(directory/'observed.npz',allow_pickle=False)
    if not np.array_equal(observations['bins_bp'],edges_bp()):
        raise ValueError('observation bins changed')
    if arrays['histograms'].shape!=(n,len(edges_bp())-1):
        raise ValueError('bank histogram dimensions')
    seeds=set(); streams=set(); families=nodes=0
    with (directory/'screen-runs.jsonl').open() as f:
        count=0
        for i,line in enumerate(f):
            r=json.loads(line); count+=1
            if i>=n or r['proposal']!=rows[i]['id'] or r['repeat']!=0:
                raise ValueError('bank log coverage mismatch')
            if r['model']!=rows[i]['model']:
                raise ValueError('bank log model mismatch')
            if not np.array_equal(r['histogram'],arrays['histograms'][i]) or not np.array_equal(r['reads'],arrays['reads'][i]):
                raise ValueError('bank log outcome mismatch')
            d=cdf_distance(arrays['reads'][i],observations['counts'])
            saved=np.inf if r['distance'] is None else r['distance']
            if d!=arrays['distances'][i] or d!=saved:
                raise ValueError('bank distance mismatch')
            seeds.add(r['seed']); streams.add(tuple(r['observation_stream']))
            families+=r['immigrants']; nodes+=r['generated_nodes']
    if count!=n or len(seeds)!=n or len(streams)!=n:
        raise ValueError('incomplete or repeated population/read streams')
    if families!=done['immigrants'] or nodes!=done['generated_nodes']:
        raise ValueError('bank work totals mismatch')
    return plan,manifest,done,arrays,observations['counts'].copy()


def observe_profile(histograms, rows, observed, seed, profile_index):
    """New independent read randomness conditional on each shared population."""
    y=np.array([observe_counts(h,int(observed.sum()),rng_for(seed,r['id'],800+profile_index))
                for h,r in zip(histograms,rows)])
    return y,np.array([cdf_distance(read,observed) for read in y])


def alpha_summary(rows, distances, immigration, threshold):
    """Equal sign-stratum allocation corresponds to alpha~Uniform(-1.5,1.5)."""
    plus=[i for i,r in enumerate(rows) if r['case']=='positive_'+immigration]
    minus=[i for i,r in enumerate(rows) if r['case']=='negative_'+immigration]
    if not plus or len(plus)!=len(minus):
        raise ValueError('equal nonempty prior allocation required for pooled alpha summaries')
    accepted=[i for i in plus+minus if np.isfinite(distances[i]) and distances[i]<=threshold]
    a=np.array([rows[i]['theta']['alpha'] for i in accepted])
    kminus=int(np.sum(a<0))
    return dict(immigration=immigration,epsilon=threshold,accepted=len(a),
                negative_accepted=kminus,positive_accepted=int(np.sum(a>=0)),
                negative_fraction=(kminus/len(a) if len(a) else None),
                alpha_quantiles025_50_975=(np.quantile(a,[.025,.5,.975]).tolist() if len(a) else None),
                accepted_proposal_ids=[rows[i]['id'] for i in accepted],
                warning=('Sparse ABC sample: fewer than 20 acceptances' if len(a)<20 else None))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bank',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--data',type=Path,default=Path(__file__).resolve().parents[2]/'data_preprocessing/mtDNA_modeling_data.csv')
    p.add_argument('--seed',type=int,default=2026101402)
    args=p.parse_args()
    if args.output.exists():
        p.error('new output directory required; never overwrite or retry silently')
    start=perf_counter()
    plan,manifest,done,bank,old_observed=load_bank(args.bank,args.data)
    if args.seed==plan['seed']:
        p.error('new read seed must differ from original bank seed')
    labels,observed=read_data(args.data,edges_bp())
    if plan['target'] not in {r['label'] for r in labels}:
        p.error('original target missing')
    original=next(i for i,r in enumerate(labels) if r['label']==plan['target'])
    if not np.array_equal(observed[original],old_observed):
        p.error('original target changed')
    args.output.mkdir(parents=True,exist_ok=False)
    source=source_info()
    source['source_sha256']['forest/examples/cross_profile_screen.py']=sha(Path(__file__))
    specification=dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),provenance=source,
        bank_directory=str(args.bank.resolve()),bank_source_commit=manifest['provenance']['git_commit'],
        bank_sha256={name:sha(args.bank/name) for name in ('plan.json','manifest.json','screen.npz','screen-runs.jsonl','screen-complete.json')},
        data_sha256=sha(args.data),seed=args.seed,original_target=plan['target'],profiles=labels,
        thresholds=list(THRESHOLDS),primary_display_threshold=.10,prior=plan['prior'],
        biological_cutoff_bp=1,observable_cutoff_bp=39.5,
        observation_stream='[seed, original proposal ID, 800 + zero-based stable profile index]',
        original_target_reads='Preserved from original bank; no replacement of its screening result.',
        interpretation='Separate marginal ABC screens conditional on immigration law. Shared population bank induces cross-profile Monte Carlo dependence. Do not multiply acceptance masks, infer a joint posterior, or equate interval overlap with universality.',
        resource_accounting='Only new multinomial read draws; no new latent populations, families or nodes.')
    dump(args.output/'plan.json',specification)
    try:
        ys=[]; ds=[]; cases=[]; alphas=[]; candidates=[]
        rows=plan['rows']
        for j,label in enumerate(labels):
            if j==original:
                y,d=bank['reads'].copy(),bank['distances'].copy()
            else:
                y,d=observe_profile(bank['histograms'],rows,observed[j],args.seed,j)
            ys.append(y);ds.append(d)
            for case in dict.fromkeys(r['case'] for r in rows):
                ix=[i for i,r in enumerate(rows) if r['case']==case]
                cases.append(dict(profile=label['label'],depth=label['n'],case=case,**success_summary(d[ix])))
            for immigration in ('uniform','beta'):
                for threshold in THRESHOLDS:
                    alphas.append(dict(profile=label['label'],depth=label['n'],**alpha_summary(rows,d,immigration,threshold)))
            candidates.append(dict(profile=label['label'],selected=select_candidates(rows,d,3)))
        np.savez_compressed(args.output/'screen.npz',reads=np.array(ys),distances=np.array(ds),observed=observed,
            bins_bp=edges_bp(),proposal_ids=np.array([r['id'] for r in rows]),profiles=np.array([r['label'] for r in labels]))
        dump(args.output/'case-summary.json',cases)
        dump(args.output/'alpha-summary.json',alphas)
        dump(args.output/'selected.json',candidates)
        complete=dict(wall_seconds=perf_counter()-start,profiles=len(labels),new_profiles=len(labels)-1,
            bank_populations=len(rows),new_populations=0,new_families=0,new_generated_nodes=0,
            new_read_draws=(len(labels)-1)*len(rows),all_profile_outcomes=len(labels)*len(rows),
            shared_bank_history_union_bound=done['history_union_bound'],
            caveat='History bound applies to shared forward bank, not the conditional ABC posterior.',
            artifacts_sha256={n:sha(args.output/n) for n in ('plan.json','screen.npz','case-summary.json','alpha-summary.json','selected.json')})
        dump(args.output/'screen-complete.json',complete)
        print(json.dumps(complete),flush=True)
    except BaseException as exc:
        dump(args.output/'FAILED.json',dict(error=repr(exc),note='Partial output diagnostic only; no outcomes filtered or replaced.'))
        raise


if __name__=='__main__':
    main()
