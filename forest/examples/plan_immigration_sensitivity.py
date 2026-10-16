"""Freeze an S04 sensitivity prior before any random population is simulated.

The original rate prior remains available. The separately authorized abundance
prior maps expected observable count to total immigration rate deterministically.
No proposal is removed based on its certificate or eventual outcome.
"""
import argparse
from dataclasses import asdict
import hashlib
from pathlib import Path
import numpy as np
from frime_forest import Immigration
from frime_forest.trees import rng_for
from abc_pilot import PRIOR, dump, model_for


def make_plan(n=1000, seed=2026101202, immigration_prior='rate',
              max_families=100000, max_total_nodes=10000000):
    from dataclasses import replace
    if immigration_prior not in ('rate','abundance'):
        raise ValueError('unknown immigration prior')
    rows=[]
    for arm,(sign,immigration) in enumerate((('positive','uniform'),('positive','beta'),
                                            ('negative','uniform'),('negative','beta'))):
        rng=rng_for(seed,arm,700)
        alpha=rng.uniform(0,1.5,n)*(1 if sign=='positive' else -1)
        boundary=np.exp(rng.uniform(np.log(50),np.log(6000),n))
        exit_rate=np.exp(rng.uniform(np.log(.5),np.log(4),n))
        bounds=(10,100) if immigration_prior=='rate' else (100,10000)
        intensity=np.exp(rng.uniform(np.log(bounds[0]),np.log(bounds[1]),n))
        shapes=np.exp(rng.uniform(np.log(.5),np.log(8),(n,2)))
        if immigration=='uniform':
            shapes[:]=1
        for i in range(n):
            theta=dict(alpha=float(alpha[i]),boundary_bp=float(boundary[i]),
                       exit_at_half_boundary=float(exit_rate[i]),immigration=float(intensity[i]),
                       immigration_a=float(shapes[i,0]),immigration_b=float(shapes[i,1]))
            model=replace(model_for([alpha[i],boundary[i],exit_rate[i],intensity[i]]),
                immigration=Immigration(immigration,a=float(shapes[i,0]),b=float(shapes[i,1])))
            mapping={}
            if immigration_prior=='abundance':
                from frime_forest.abundance_prior import with_expected_abundance
                model,coefficient=with_expected_abundance(model,float(intensity[i]))
                theta['expected_observable_count']=float(intensity[i])
                theta['immigration']=model.c_i
                mapping['abundance_mapping']=asdict(coefficient)
            rows.append(dict(id=arm*n+i,case=sign+'_'+immigration,theta=theta,
                             model=asdict(model),**mapping))
            if (i+1)%100==0:
                print(f'{sign}_{immigration}: {i+1}/{n} prior settings mapped',flush=True)
    prior={k:v for k,v in PRIOR.items() if k!='alpha'}
    if immigration_prior=='abundance':
        prior.pop('immigration',None)
        prior['expected_observable_count']=['log_uniform',100,10000]
    return dict(seed=seed,target='S04',delta=1e-10,top_per_case=3,confirmation_repeats=100,
        max_families=max_families,max_nodes=1000000,max_total_nodes=max_total_nodes,
        immigration_prior=immigration_prior,
        prior=dict(alpha_negative=['uniform',-1.5,0],alpha_positive=['uniform',0,1.5],
            immigration_a=['log_uniform',.5,8],immigration_b=['log_uniform',.5,8],
            **prior),
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        notes=['Four separate independent parameter banks; no case posterior probabilities.',
               ('Immigration coordinate is total arrivals per boundary fragmentation time, rho(B)=1.'
                if immigration_prior=='rate' else
                'Approved 2026-10-14: shared prior on expected observable abundance; C_I=mu/K. '
                'K is numerically integrated and refinement checked; not a history certificate. '
                'Abundance is latent population size, not read depth or measured concentration.'),
               'Preflight failure stops the entire plan before simulations; no expensive proposals filtered.',
               'Top three finite-distance candidates per case selected by screen CDF gap, then proposal ID.',
               'Common descriptive distance thresholds .05,.10,.15; fresh 100-population candidate repeats.',
               'Original pilot is not a reproduction of the thesis fitted regimes.'],rows=rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--proposals',type=int,default=1000)
    p.add_argument('--seed',type=int,default=2026101202)
    p.add_argument('--immigration-prior',choices=('rate','abundance'),default='rate')
    p.add_argument('--max-families',type=int,default=100000)
    p.add_argument('--max-total-nodes',type=int,default=10000000)
    args=p.parse_args()
    if args.output.exists() or args.proposals<1:
        p.error('positive proposal count and a new output path required')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    dump(args.output,make_plan(args.proposals,args.seed,args.immigration_prior,
                              args.max_families,args.max_total_nodes))


if __name__=='__main__':
    main()
