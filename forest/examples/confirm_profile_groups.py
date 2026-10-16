"""Fresh fixed-candidate confirmations for the six selected mtDNA profiles.

Prepare freezes candidates and certifies resources; run generates independent
populations per profile/candidate/repeat. No refitting or common-alpha test here.
"""
from __future__ import annotations
import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import json
import os
from pathlib import Path
from time import perf_counter
import numpy as np
from frime_forest import MonteCarlo, replicate_seed
from frime_forest.abc import cdf_distance, observe_counts
from frime_forest.certificate import certify
from frime_forest.trees import rng_for
from abc_pilot import dump, edges_bp, read_data, L_BP
from cross_profile_screen import load_bank, sha
from immigration_sensitivity import model_from, source_info, select_candidates, predictive_summary

PROFILES=('S01','S05','S06','S07','S08','S09')
CASES=('positive_beta','negative_beta')
LATENT_COMPONENT=900
READ_COMPONENT=901


def sources():
    info=source_info()
    for name in ('cross_profile_screen.py','confirm_profile_groups.py'):
        p=Path(__file__).with_name(name)
        info['source_sha256']['forest/examples/'+name]=sha(p)
    return info


def select_cells(rows, distances, labels, profiles=PROFILES, top=3):
    cells=[]
    for profile in profiles:
        j=labels.index(profile)
        selected=select_candidates(rows,distances[j],top)
        for case in CASES:
            case_rows=[r for r in selected if r['case']==case]
            if len(case_rows)!=top:
                raise ValueError('insufficient finite screened candidates; do not substitute')
            for rank,row in enumerate(case_rows,1):
                cells.append(dict(cell=len(cells),profile=profile,rank=rank,**row))
    return cells


def stream_index(cell, repeat, repeats):
    if type(cell) is not int or cell<0 or type(repeat) is not int or not 0<=repeat<repeats:
        raise ValueError('invalid stream index')
    return cell*repeats+repeat


def prepare(screen, output, data, seed, repeats):
    if type(repeats) is not int or repeats<2:
        raise ValueError('at least two independent repeats required')
    spec=json.loads((screen/'plan.json').read_text())
    done=json.loads((screen/'screen-complete.json').read_text())
    for name,h in done['artifacts_sha256'].items():
        if sha(screen/name)!=h:raise ValueError('cross-profile screen artifact changed')
    bank=Path(spec['bank_directory'])
    original,manifest,_,_,_=load_bank(bank,data)
    for name,h in spec['bank_sha256'].items():
        if sha(bank/name)!=h:raise ValueError('original bank changed')
    if seed in (original['seed'],spec['seed']):raise ValueError('new root seed required')
    labels,observed=read_data(data,edges_bp()); names=[r['label'] for r in labels]
    with np.load(screen/'screen.npz',allow_pickle=False) as z:
        if names!=z['profiles'].tolist() or not np.array_equal(observed,z['observed']):
            raise ValueError('observations changed')
        cells=select_cells(original['rows'],z['distances'],names)
    saved=json.loads((screen/'selected.json').read_text())
    for c in cells:
        row={k:v for k,v in c.items() if k not in ('cell','profile','rank')}
        if row not in next(x['selected'] for x in saved if x['profile']==c['profile']):
            raise ValueError('candidate selection differs from saved screening rule')
    certificates=[]
    for c in cells:
        model=model_from(c['model']);cert=certify(model,original['delta'])
        certificates.append(dict(cell=c['cell'],expected_families=model.c_i*cert.T,certificate=asdict(cert)))
    if any(c['expected_families']>original['max_families'] for c in certificates):
        raise ValueError('entire confirmation plan exceeds cap; no rows filtered')
    plan=dict(seed=seed,repeats=repeats,profiles=list(PROFILES),cells=cells,
        delta=original['delta'],max_families=original['max_families'],max_nodes=original['max_nodes'],
        max_total_nodes=original['max_total_nodes'],certificates=certificates,
        expected_total_families=repeats*sum(c['expected_families'] for c in certificates),
        profile_metadata=[r for r in labels if r['label'] in PROFILES],
        observed={r['label']:observed[j].tolist() for j,r in enumerate(labels) if r['label'] in PROFILES},
        screen_directory=str(screen.resolve()),bank_directory=str(bank.resolve()),
        data_sha256=sha(data),input_sha256={str((screen/name).resolve()):sha(screen/name) for name in ('plan.json','screen.npz','selected.json','screen-complete.json')},
        provenance=sources(),timestamp_utc=datetime.now(timezone.utc).isoformat(),
        observation_stream='[root seed, cell*repeats+repeat, 901]',
        population_stream='replicate_seed(replicate_seed(root seed,900),cell*repeats+repeat)',
        primary_threshold=.10,thresholds=[.05,.10,.15],
        analysis=dict(read_bootstrap_repeats=1000,read_bootstrap_component=902,
            primary_candidates='Screen rank 1 per sign/profile; ranks 2/3 are prespecified robustness checks.',
            recovery='Refit every fresh readout using the full fixed 2000-proposal Beta bank at its profile depth; retain no-acceptance cases.',
            bootstrap='Resample observed multinomial reads at original depth, refit same fixed bank. Conditional read-sampling sensitivity, not calibrated biological confidence intervals.',
            variance='Separate latent population CDF variance from expected multinomial read variance by total variance.',
            caveat='Selected fixed-parameter checks; not a shared-alpha test or global calibration. No new candidate selected as primary after confirmation.'))
    output.mkdir(parents=True,exist_ok=False)
    dump(output/'plan.json',plan)
    dump(output/'prepare-complete.json',dict(plan_sha256=sha(output/'plan.json'),cells=len(cells),
        populations=len(cells)*repeats,expected_total_families=plan['expected_total_families']))
    return plan


def check_plan(output, data):
    ready=json.loads((output/'prepare-complete.json').read_text())
    if sha(output/'plan.json')!=ready['plan_sha256']:raise ValueError('frozen plan changed')
    plan=json.loads((output/'plan.json').read_text())
    if sha(data)!=plan['data_sha256']:raise ValueError('clinical data changed')
    for name,h in plan['input_sha256'].items():
        if sha(Path(name))!=h:raise ValueError('screen input changed')
    if sources()['source_sha256']!=plan['provenance']['source_sha256']:
        raise ValueError('scientific source changed')
    return plan


def run(plan, output, workers, batch_size):
    started=output/'run-started.json'
    if started.exists():raise ValueError('run already attempted; never silently retry or overwrite')
    dump(started,dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),pid=os.getpid(),
        provenance=sources(),workers=workers,batch_size=batch_size))
    start=perf_counter();total=dict(populations=0,immigrants=0,generated_nodes=0,history_union_bound=0.)
    hs=[];ys=[];indices=[];summaries=[];distances=[]
    ensemble=replicate_seed(plan['seed'],LATENT_COMPONENT)
    try:
        with (output/'runs.jsonl').open('x') as log,MonteCarlo(workers=workers,batch_size=batch_size,
            max_families=plan['max_families'],max_nodes=plan['max_nodes'],max_total_nodes=plan['max_total_nodes']) as mc:
            for cell in plan['cells']:
                local_y=[];observed=np.array(plan['observed'][cell['profile']]);model=model_from(cell['model'])
                for repeat in range(plan['repeats']):
                    index=stream_index(cell['cell'],repeat,plan['repeats'])
                    result=mc.run(model,1,bins=edges_bp()/L_BP,seed=ensemble,start_replicate=index,delta=plan['delta'])
                    h=result.histograms[0]
                    y=observe_counts(h,int(observed.sum()),rng_for(plan['seed'],index,READ_COMPONENT))
                    d=cdf_distance(y,observed)
                    rec=dict(cell=cell['cell'],profile=cell['profile'],case=cell['case'],rank=cell['rank'],proposal=cell['id'],repeat=repeat,
                        theta=cell['theta'],model=cell['model'],certificate=asdict(result.certificate),
                        started_utc=result.started_utc,seed=str(result.seeds[0]),ensemble_seed=str(ensemble),
                        observation_stream=[plan['seed'],index,READ_COMPONENT],workers=workers,batch_size=batch_size,mode='lazy',
                        histogram=h.tolist(),reads=y.tolist(),distance=d if np.isfinite(d) else None,
                        immigrants=int(result.immigrants[0]),generated_nodes=int(result.generated_nodes[0]),
                        observable_count=int(result.counts[0]),observable_mass=float(result.masses[0]),seconds=float(result.seconds[0]))
                    log.write(json.dumps(rec,allow_nan=False)+'\n');log.flush()
                    hs.append(h);ys.append(y);local_y.append(y);indices.append((cell['cell'],repeat));distances.append(d)
                    total['populations']+=1;total['immigrants']+=rec['immigrants'];total['generated_nodes']+=rec['generated_nodes']
                    total['history_union_bound']+=result.certificate.lambda_bound
                    if total['populations']%10==0:
                        dump(output/'progress.json',dict(**total,profile=cell['profile'],cell=cell['cell'],
                            wall_seconds=perf_counter()-start,total_planned=len(plan['cells'])*plan['repeats']))
                s=predictive_summary(np.array(local_y),observed)
                s.update(cell=cell['cell'],profile=cell['profile'],case=cell['case'],rank=cell['rank'],proposal=cell['id'],theta=cell['theta'])
                summaries.append(s);dump(output/'partial-summary.json',summaries)
                print(json.dumps(dict(profile=cell['profile'],case=cell['case'],rank=cell['rank'],
                    completed=total['populations'],total=len(plan['cells'])*plan['repeats'],wall_seconds=perf_counter()-start)),flush=True)
        np.savez_compressed(output/'confirm.npz',histograms=np.array(hs),reads=np.array(ys),indices=np.array(indices),distances=np.array(distances))
        dump(output/'summary.json',summaries)
        total.update(wall_seconds=perf_counter()-start,finished_utc=datetime.now(timezone.utc).isoformat(),
            artifacts_sha256={name:sha(output/name) for name in ('plan.json','runs.jsonl','confirm.npz','summary.json')})
        dump(output/'complete.json',total);print(json.dumps(total),flush=True)
        return total
    except BaseException as exc:
        dump(output/'FAILED.json',dict(error=repr(exc),**total,wall_seconds=perf_counter()-start,
            note='Partial diagnostics only. No failed draws or expensive proposals replaced.'))
        raise


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--stage',choices=('prepare','run'),required=True)
    p.add_argument('--screen',type=Path)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--data',type=Path,default=Path(__file__).resolve().parents[2]/'data_preprocessing/mtDNA_modeling_data.csv')
    p.add_argument('--seed',type=int,default=2026101601)
    p.add_argument('--repeats',type=int,default=100)
    p.add_argument('--workers',type=int,default=8)
    p.add_argument('--batch-size',type=int,default=128)
    args=p.parse_args()
    if args.stage=='prepare':
        if args.screen is None:p.error('--screen required')
        plan=prepare(args.screen,args.output,args.data,args.seed,args.repeats)
        print(json.dumps(dict(cells=len(plan['cells']),populations=len(plan['cells'])*plan['repeats'],expected_total_families=plan['expected_total_families'])))
    else:run(check_plan(args.output,args.data),args.output,args.workers,args.batch_size)


if __name__=='__main__':main()
