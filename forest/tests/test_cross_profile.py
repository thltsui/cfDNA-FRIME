import importlib
from pathlib import Path
import numpy as np
import pytest


@pytest.fixture
def cross(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1]/'examples'))
    return importlib.import_module('cross_profile_screen')


def test_profile_read_law_preserves_failure_and_depth(cross):
    rows=[{'id':i} for i in range(3001)]
    h=np.tile([1,3],(len(rows),1));h[0]=0
    y,d=cross.observe_profile(h,rows,np.array([5,15]),21,2)
    assert y[0].sum()==0 and np.isinf(d[0])
    assert np.all(y[1:].sum(axis=1)==20)
    # Analytic mean and MCSE for independent Binomial(20,.25) read counts.
    assert abs(y[1:,0].mean()-5)<5*np.sqrt(20*.25*.75/3000)
    other,_=cross.observe_profile(h,rows,np.array([5,15]),21,3)
    again,_=cross.observe_profile(h,rows,np.array([5,15]),21,2)
    assert np.array_equal(y,again)
    assert not np.array_equal(y,other)


def test_alpha_screen_includes_failed_draws_in_prior_bank(cross):
    rows=[dict(id=i,case=c,theta={'alpha':a}) for i,(c,a) in enumerate([
        ('positive_beta',.5),('positive_beta',1.),('negative_beta',-.5),('negative_beta',-1.)])]
    s=cross.alpha_summary(rows,np.array([.02,np.inf,.02,.2]),'beta',.05)
    assert s['accepted_proposal_ids']==[0,2]
    assert s['negative_fraction']==.5 and s['accepted']==2
    assert s['alpha_quantiles025_50_975'][1]==0
    empty=cross.alpha_summary(rows,np.full(4,np.inf),'beta',.05)
    assert empty['accepted']==0 and empty['negative_fraction'] is None
    with pytest.raises(ValueError,match='equal nonempty'):
        cross.alpha_summary(rows[:-1],np.zeros(3),'beta',.05)


def test_complete_bank_required_and_hash_guard(cross,tmp_path,monkeypatch):
    import json
    dump=lambda p,x:p.write_text(json.dumps(x))
    data=tmp_path/'data.csv';data.write_text('fixture')
    out=tmp_path/'bank';out.mkdir()
    dump(out/'plan.json',{'rows':[{'id':0,'model':{}}]})
    np.savez_compressed(out/'screen.npz',histograms=np.ones((1,56),dtype=int),reads=np.ones((1,56),dtype=int),distances=[0.],indices=[[0,0]])
    dump(out/'manifest.json',dict(plan_sha256=cross.sha(out/'plan.json'),data_sha256=cross.sha(data),observation_sha256={},provenance={'source_sha256':{}}))
    dump(out/'screen-complete.json',dict(populations=0,artifacts_sha256={'screen.npz':cross.sha(out/'screen.npz')}))
    monkeypatch.setattr(cross,'source_info',lambda:{'source_sha256':{}})
    with pytest.raises(ValueError,match='complete original'):
        cross.load_bank(out,data)
    (out/'screen.npz').write_bytes(b'changed')
    with pytest.raises(ValueError,match='artifact changed'):
        cross.load_bank(out,data)
