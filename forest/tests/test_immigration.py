"""Independent checks of immigrant marks, compatibility and stationary laws."""
from dataclasses import replace, asdict
import math
import json

import numpy as np
import pytest

from frime_forest import (Immigration, Model, Exit, certify, make_family,
    make_skeleton, make_plan, presample, assemble, MonteCarlo, ResourceLimitError)
from frime_forest.trees import rng_for


@pytest.mark.parametrize('kwargs',[dict(kind='normal'),dict(a=0),dict(b=-1),
    dict(kind='beta',a=float('inf')),dict(kind='beta',a=1e308,b=1e308),dict(a=2)])
def test_invalid_immigration(kwargs):
    with pytest.raises(ValueError):
        Immigration(**kwargs)


def test_beta11_preserves_uniform_marks_certificates_and_bank():
    a=Model(L=2,ell=.2,c_i=1)
    b=replace(a,immigration=Immigration('beta'))
    for i in range(20):
        assert a.sample_root(rng_for(9,i,2))==b.sample_root(rng_for(9,i,2))
    ca,cb=certify(a),certify(b)
    assert (ca.T,ca.log_A,ca.kappa,ca.lambda_bound)==(cb.T,cb.log_A,cb.kappa,cb.lambda_bound)
    plan=make_plan(b,seed=23)
    bank=presample(a,len(plan.arrivals),seed=23)
    assert len(assemble(plan,bank).families)==len(bank)


@pytest.mark.parametrize('a,b',[(2,1),(1,2),(.5,4),(4,.5)])
def test_root_moments_without_cutoff_conditioning(a,b):
    model=Model(L=3,ell=2,immigration=Immigration('beta',a,b))
    values=np.array([model.sample_root(rng_for(4612,i,2))/3 for i in range(5000)])
    for p in (1,2):
        exact=math.prod((a+j)/(a+b+j) for j in range(p))
        samples=values**p
        z=(samples.mean()-exact)/(samples.std(ddof=1)/math.sqrt(len(samples)))
        assert abs(z)<6, (a,b,p,z)
    assert np.any(values<model.ell/model.L)
    for i in range(30):
        family=make_family(model,i,4612,0)
        assert bool(len(family.sizes))==(family.root_size>model.ell)


@pytest.mark.parametrize('q',[0.,1.,float('nan'),float('inf')])
def test_beta_endpoint_raises_without_retry(q):
    class Fixed:
        calls=0
        def beta(self,*args):
            self.calls+=1
            return q
    rng=Fixed()
    with pytest.raises(FloatingPointError):
        Model(immigration=Immigration('beta',2,1)).sample_root(rng)
    assert rng.calls==1


@pytest.mark.parametrize('a,b',[(2,1),(1,2),(.5,8),(8,.5)])
def test_full_beta_moment_dominates_truncated_moment(a,b):
    model=Model(L=3,ell=.4,immigration=Immigration('beta',a,b))
    for p in (2,3,4,6,8):
        c=certify(model,1e-9,powers=(p,))
        exact_log=p*math.log(3/.4)+math.log(math.prod((a+j)/(a+b+j) for j in range(p)))
        assert exact_log<=c.log_A<=exact_log+1e-8
        assert c.prefactor_method=='beta_full_moment_upper_bound'
        assert 0<c.lambda_bound<=1e-9*(1+1e-12)


def test_bank_cannot_change_root_distribution():
    model=Model(c_i=1,ell=.1,immigration=Immigration('beta',2,1))
    other=replace(model,immigration=Immigration('beta',1,2))
    tree=make_skeleton(model,0,91)
    with pytest.raises(ValueError,match='immigrant law'):
        tree.prune(other)
    plan=make_plan(other,seed=22)
    assert len(plan.arrivals)>0
    bank=presample(model,len(plan.arrivals),seed=22)
    with pytest.raises(ValueError,match='immigrant law'):
        assemble(plan,bank)


@pytest.mark.parametrize('mode',['lazy','skeleton'])
def test_beta_population_worker_batch_reproducibility(mode):
    model=Model(ell=.1,c_i=2,alpha=-.5,immigration=Immigration('beta',.5,3))
    kw=dict(bins=[.1,.2,.5,1],seed=22333,delta=1e-7,mode=mode)
    with MonteCarlo(workers=1,batch_size=1) as mc:
        serial=mc.run(model,4,**kw)
    with MonteCarlo(workers=2,batch_size=7) as mc:
        parallel=mc.run(model,4,**kw)
    for name in ('histograms','counts','masses','immigrants','generated_nodes'):
        np.testing.assert_array_equal(getattr(serial,name),getattr(parallel,name))
    with MonteCarlo(workers=1,max_families=1) as mc:
        with pytest.raises(ResourceLimitError):
            mc.run(model,1,**kw)


@pytest.mark.parametrize('mode',['lazy','skeleton'])
@pytest.mark.parametrize('a,b,sign',[(2,1,1),(1,2,-1)])
def test_beta_stationary_first_moment_with_exit(mode,a,b,sign,record_property):
    model=Model(L=2,ell=.2,c_f=1,alpha=0,c_i=4,exit=Exit(c_e=1),
                immigration=Immigration('beta',a,b))
    bins=np.array([.2,.5,1,2])
    expected=2*np.log(bins[1:]/bins[:-1])+sign/4*(bins[1:]**2-bins[:-1]**2)
    with MonteCarlo() as mc:
        r=mc.run(model,300,bins=bins,seed=65231+(10 if sign<0 else 0)+(1 if mode=='skeleton' else 0),
                 delta=1e-10,mode=mode)
    mean=r.histograms.mean(axis=0)
    se=r.histograms.std(axis=0,ddof=1)/math.sqrt(300)
    z=(mean-expected)/se
    record_property("statistical_validation",json.dumps(dict(model=asdict(model),seed=r.seed,
        n=300,mode=r.mode,certificate=asdict(r.certificate),mean=mean.tolist(),
        expected=expected.tolist(),mcse=se.tolist(),z=z.tolist(),
        generated_nodes=int(r.generated_nodes.sum()),immigrants=int(r.immigrants.sum()),
        runtime_seconds=r.elapsed_seconds)))
    assert np.all(np.abs(z)<6),dict(mean=mean,analytic=expected,mcse=se,z=z)


@pytest.mark.parametrize('a,b,coefficient,sign',[(2,1,4/3,1),(1,2,2/3,-1)])
def test_negative_alpha_beta_first_moment(a,b,coefficient,sign,record_property):
    model=Model(L=2,ell=.2,c_f=2,alpha=-1,c_i=2,exit=Exit(c_e=0),
                immigration=Immigration('beta',a,b))
    bins=np.array([.2,.5,1,2])
    expected=coefficient*2*np.log(bins[1:]/bins[:-1])+sign/18*(bins[1:]**3-bins[:-1]**3)
    with MonteCarlo() as mc:
        r=mc.run(model,300,bins=bins,seed=736411+a,delta=1e-10)
    mean=r.histograms.mean(axis=0);se=r.histograms.std(axis=0,ddof=1)/math.sqrt(300)
    z=(mean-expected)/se
    record_property("statistical_validation",json.dumps(dict(model=asdict(model),seed=r.seed,
        n=300,mode=r.mode,certificate=asdict(r.certificate),mean=mean.tolist(),
        expected=expected.tolist(),mcse=se.tolist(),z=z.tolist(),
        generated_nodes=int(r.generated_nodes.sum()),immigrants=int(r.immigrants.sum()),
        runtime_seconds=r.elapsed_seconds)))
    assert np.all(np.abs(z)<6),dict(mean=mean,analytic=expected,mcse=se,z=z)


def test_numpy_integer_shapes_are_normalized_before_sum():
    law=Immigration("beta",np.int64(2**62),np.int64(2**62))
    assert type(law.a) is float and type(law.b) is float
    assert law.a+law.b==float(2**63)
