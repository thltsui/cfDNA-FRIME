"""Freeze the originally proposed S04 sensitivity prior before any simulation.

This rate prior is retained for an honest resource audit even if it proves too
expensive. No proposal is removed based on its certificate or eventual outcome.
"""
import argparse
from dataclasses import asdict
import hashlib
from pathlib import Path
import numpy as np
from frime_forest import Immigration
from frime_forest.trees import rng_for
from abc_pilot import PRIOR, dump, model_for


def make_plan(n=1000, seed=2026101202):
    from dataclasses import replace
    rows=[]
    for arm,(sign,immigration) in enumerate((('positive','uniform'),('positive','beta'),
                                            ('negative','uniform'),('negative','beta'))):
        rng=rng_for(seed,arm,700)
        alpha=rng.uniform(0,1.5,n)*(1 if sign=='positive' else -1)
        boundary=np.exp(rng.uniform(np.log(50),np.log(6000),n))
        exit_rate=np.exp(rng.uniform(np.log(.5),np.log(4),n))
        intensity=np.exp(rng.uniform(np.log(10),np.log(100),n))
        shapes=np.exp(rng.uniform(np.log(.5),np.log(8),(n,2)))
        if immigration=='uniform':
            shapes[:]=1
        for i in range(n):
            theta=dict(alpha=float(alpha[i]),boundary_bp=float(boundary[i]),
                       exit_at_half_boundary=float(exit_rate[i]),immigration=float(intensity[i]),
                       immigration_a=float(shapes[i,0]),immigration_b=float(shapes[i,1]))
            model=replace(model_for([alpha[i],boundary[i],exit_rate[i],intensity[i]]),
                immigration=Immigration(immigration,a=float(shapes[i,0]),b=float(shapes[i,1])))
            rows.append(dict(id=arm*n+i,case=sign+'_'+immigration,theta=theta,model=asdict(model)))
    return dict(seed=seed,target='S04',delta=1e-10,top_per_case=3,confirmation_repeats=100,
        max_families=100000,max_nodes=1000000,max_total_nodes=10000000,
        prior=dict(alpha_negative=['uniform',-1.5,0],alpha_positive=['uniform',0,1.5],
            immigration_a=['log_uniform',.5,8],immigration_b=['log_uniform',.5,8],
            **{k:v for k,v in PRIOR.items() if k!='alpha'}),
        generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        notes=['Four separate independent parameter banks; no case posterior probabilities.',
               'Immigration coordinate is total arrivals per boundary fragmentation time, rho(B)=1.',
               'Preflight failure stops the entire plan before simulations; no expensive proposals filtered.',
               'Top three finite-distance candidates per case selected by screen CDF gap, then proposal ID.',
               'Common descriptive distance thresholds .05,.10,.15; fresh 100-population candidate repeats.',
               'Original pilot is not a reproduction of the thesis fitted regimes.'],rows=rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--proposals',type=int,default=1000)
    p.add_argument('--seed',type=int,default=2026101202)
    args=p.parse_args()
    if args.output.exists() or args.proposals<1:
        p.error('positive proposal count and a new output path required')
    args.output.parent.mkdir(parents=True,exist_ok=True)
    dump(args.output,make_plan(args.proposals,args.seed))


if __name__=='__main__':
    main()
