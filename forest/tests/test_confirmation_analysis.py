import importlib.util
from pathlib import Path
import numpy as np

spec=importlib.util.spec_from_file_location('analysis',Path(__file__).resolve().parents[1]/'examples'/'analyze_profile_confirmation.py')
a=importlib.util.module_from_spec(spec);spec.loader.exec_module(a)

def test_refit_matches_scalar_rejection_and_retains_failures():
    rng=np.random.default_rng(31)
    bank=rng.multinomial(30,[.2,.3,.5],size=50);bank[0]=0
    targets=rng.multinomial(30,[.2,.3,.5],size=19);targets[0]=0
    alpha=np.linspace(-1.5,1.5,len(bank));got=a.refit(bank,targets,alpha)
    for i,t in enumerate(targets):
        d=np.array([np.max(np.abs(b.cumsum()/b.sum()-t.cumsum()/t.sum())) if b.sum() and t.sum() else np.inf for b in bank])
        for j,e in enumerate(a.EPS):
            values=alpha[d<=e]
            assert got['accepted'][i,j]==len(values)
            if len(values):
                assert got['median'][i,j]==np.median(values)
                assert got['negative_fraction'][i,j]==np.mean(values<0)
            else:assert np.isnan(got['median'][i,j])

def test_total_variance_decomposition_known_mixture():
    rng=np.random.default_rng(32);n=100
    # Two equiprobable latent probabilities .2 and .8: variance .09;
    # mean conditional read variance .16/n.
    h=np.tile([[2,8],[8,2]],(10000,1))
    y=np.array([rng.multinomial(n,r/r.sum()) for r in h])
    v=a.variance_components(h,y,n)
    assert np.isclose(v['latent'][0],.09*20000/19999)
    assert np.isclose(v['sampling'][0],.0016)
    assert abs(v['empirical'][0]-v['total'][0])<.002
    assert v['total'][-1]==0

def test_zero_read_target_is_not_inferred_from_zero_bank():
    r=a.refit(np.array([[0,0],[10,0]]),np.array([[0,0]]),np.array([-1.,1.]))
    assert not r['accepted'].any()
    assert np.isnan(r['median']).all()
