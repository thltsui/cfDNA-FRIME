import importlib
from dataclasses import asdict
import json
from pathlib import Path
import numpy as np
import pytest
from frime_forest import Model,Exit,replicate_seed
from frime_forest.trees import ResourceLimitError

@pytest.fixture
def study(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1]/'examples'))
    return importlib.import_module('confirm_profile_groups')

def test_profile_selection_is_ranked_before_confirmation(study):
    rows=[dict(id=i,case=c) for i,c in enumerate(['positive_beta']*4+['negative_beta']*4+['positive_uniform']*4)]
    cells=study.select_cells(rows,np.array([[.3,.1,.1,.2,.3,.2,np.inf,.1,0,0,0,0]]),['S01'],('S01',),3)
    assert [r['id'] for r in cells]==[1,2,3,7,5,4]
    assert [r['rank'] for r in cells]==[1,2,3,1,2,3]
    assert [r['cell'] for r in cells]==list(range(6))
    with pytest.raises(ValueError,match='insufficient'):
        study.select_cells(rows,np.full((1,12),np.inf),['S01'],('S01',),3)

def test_streams_distinct_even_for_same_candidate_across_profiles(study):
    ids=[study.stream_index(c,j,100) for c in range(36) for j in range(100)]
    assert len(set(ids))==3600
    root=replicate_seed(2026101601,900)
    new={replicate_seed(root,i) for i in ids}
    old={replicate_seed(replicate_seed(2026101401,720),i) for i in ids}
    assert not new & old
    with pytest.raises(ValueError):study.stream_index(0,100,100)

def test_real_fresh_run_logs_complete_populations_and_refuses_retry(study,tmp_path):
    model=asdict(Model(L=1,ell=39.5/16569,c_f=1,alpha=0,c_i=2,exit=Exit(c_e=2)))
    plan=dict(seed=44,repeats=3,delta=1e-7,max_families=10000,max_nodes=10000,max_total_nodes=100000,
        observed={'S01':[20]+[0]*55,'S05':[30]+[0]*55},
        cells=[dict(cell=i,profile=p,case='positive_beta',rank=1,id=55,theta={},model=model) for i,p in enumerate(['S01','S05'])])
    (tmp_path/'plan.json').write_text(json.dumps(plan))
    result=study.run(plan,tmp_path,1,16)
    logs=[json.loads(s) for s in (tmp_path/'runs.jsonl').read_text().splitlines()]
    assert result['populations']==6
    assert len({r['seed'] for r in logs})==6
    assert len({tuple(r['observation_stream']) for r in logs})==6
    assert all(sum(r['reads']) in (0,20 if r['profile']=='S01' else 30) for r in logs)
    with pytest.raises(ValueError,match='already attempted'):study.run(plan,tmp_path,1,16)

def test_resource_failure_is_not_a_completed_partial_run(study,tmp_path):
    model=asdict(Model(L=1,ell=.01,c_f=1,alpha=0,c_i=100,exit=Exit(c_e=1)))
    plan=dict(seed=45,repeats=2,delta=1e-7,max_families=1,max_nodes=10000,max_total_nodes=100000,
        observed={'S01':[20]+[0]*55},cells=[dict(cell=0,profile='S01',case='positive_beta',rank=1,id=1,theta={},model=model)])
    with pytest.raises(ResourceLimitError):study.run(plan,tmp_path,1,16)
    assert (tmp_path/'FAILED.json').exists()
    assert not (tmp_path/'complete.json').exists()
