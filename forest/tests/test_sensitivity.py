"""Frozen-plan integrity and fresh-simulation checks for the study executor."""
import importlib
import json
from pathlib import Path
import sys

import numpy as np
import pytest


@pytest.fixture
def study(monkeypatch):
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1]/'examples'))
    return importlib.import_module('immigration_sensitivity')


@pytest.fixture
def plan_files(tmp_path):
    from dataclasses import asdict
    from frime_forest import Model, Exit
    model=asdict(Model(L=1,ell=39.5/16569,c_f=1,alpha=0,c_i=5,exit=Exit(c_e=1)))
    plan=dict(seed=2026101203,target='S04',delta=1e-8,top_per_case=1,
              confirmation_repeats=2,max_families=100000,max_nodes=1000000,
              max_total_nodes=10000000,
              rows=[dict(id=i,case=f'case{i}',theta={},model=model) for i in range(4)])
    data=tmp_path/'data.csv'
    data.write_text('id,chr,length,freq\n'+''.join(f'{i},chrM,100,50\n' for i in range(4)))
    path=tmp_path/'plan.json';path.write_text(json.dumps(plan))
    return path,data,tmp_path/'out'


def invoke(study,monkeypatch,files,stage):
    plan,data,out=files
    argv=['study','--stage',stage,'--data',str(data),'--output',str(out),'--workers','1']
    if stage=='check':
        argv+=['--plan',str(plan)]
    monkeypatch.setattr(sys,'argv',argv)
    study.main()


def test_selection_ties_and_no_read_failures(study):
    rows=[dict(id=i,case=c) for i,c in ((3,'a'),(1,'a'),(4,'a'),(8,'b'),(9,'b'))]
    chosen=study.select_candidates(rows,[.1,.1,float('inf'),.2,float('inf')],2)
    assert [r['id'] for r in chosen]==[1,3,8]


@pytest.mark.parametrize('key,value',[('confirmation_repeats',0),('top_per_case',0),
                                      ('confirmation_repeats',True)])
def test_invalid_plan_cannot_reuse_stream_zero(study,monkeypatch,plan_files,key,value):
    p=plan_files[0];plan=json.loads(p.read_text());plan[key]=value;p.write_text(json.dumps(plan))
    with pytest.raises(SystemExit):
        invoke(study,monkeypatch,plan_files,'check')
    assert not (plan_files[2]/'check-complete.json').exists()


def test_fresh_confirmation_and_artifact_integrity(study,monkeypatch,plan_files):
    invoke(study,monkeypatch,plan_files,'check')
    invoke(study,monkeypatch,plan_files,'screen')
    out=plan_files[2]
    selected=out/'selected.json';saved=selected.read_bytes()
    rows=json.loads(saved);rows[0]['model']['c_i']=100
    selected.write_text(json.dumps(rows))
    with pytest.raises(SystemExit):
        invoke(study,monkeypatch,plan_files,'confirm')
    assert not (out/'confirm-started.json').exists()
    selected.write_bytes(saved)
    invoke(study,monkeypatch,plan_files,'confirm')
    screen=[json.loads(s) for s in (out/'screen-runs.jsonl').read_text().splitlines()]
    confirm=[json.loads(s) for s in (out/'confirm-runs.jsonl').read_text().splitlines()]
    seeds=[r['seed'] for r in screen+confirm]
    assert len(screen)==4 and len(confirm)==8 and len(set(seeds))==12
    observations=[tuple(r['observation_stream']) for r in screen+confirm]
    assert len(set(observations))==12
    assert all(sum(r['reads'])==50 for r in screen+confirm)
    result=json.loads((out/'confirm-complete.json').read_text())
    assert result['populations']==8


def test_observation_change_and_resource_failure_stop_whole_stage(study,monkeypatch,plan_files):
    invoke(study,monkeypatch,plan_files,'check')
    out=plan_files[2]
    np.savez_compressed(out/'observed.npz',counts=[50,0],bins_bp=[0,1,2])
    with pytest.raises(SystemExit):
        invoke(study,monkeypatch,plan_files,'screen')
    assert not (out/'screen-started.json').exists()


def test_resource_preflight_does_not_keep_affordable_subset(study,monkeypatch,plan_files):
    p=plan_files[0];plan=json.loads(p.read_text());plan['max_families']=1
    p.write_text(json.dumps(plan))
    with pytest.raises(RuntimeError,match='no proposals filtered'):
        invoke(study,monkeypatch,plan_files,'check')
    out=plan_files[2]
    assert len(json.loads((out/'certificates.json').read_text()))==4
    assert (out/'check-FAILED.json').exists()
    assert not (out/'check-complete.json').exists()
    assert not list(out.glob('*-runs.jsonl'))
