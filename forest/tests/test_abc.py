import numpy as np
import pytest
from frime_forest.abc import observe_counts, cdf_distance, rejection_indices


def test_empty_population_is_no_read_outcome_not_resampled():
    class NoRandomCalls:
        def multinomial(self, *args):
            raise AssertionError("empty population must not draw")
    h = observe_counts([0, 0, 0], 100, NoRandomCalls())
    assert np.array_equal(h, [0, 0, 0])
    assert np.isinf(cdf_distance(h, [10, 20, 70]))
    assert np.array_equal(rejection_indices([float("inf"), .1], 1), [1])


def test_observation_uses_fragment_counts_not_length_mass():
    class Capture:
        def multinomial(self, n, p):
            assert n == 100
            assert np.array_equal(p, [.25, 0, .75])
            return np.array([25, 0, 75])
    assert observe_counts([1, 0, 3], 100, Capture()).sum() == 100


def test_cdf_distance_and_inclusive_fixed_tolerance():
    assert cdf_distance([1, 1, 2], [0, 2, 2]) == .25
    assert cdf_distance([1, 0], [0, 1]) == 1
    assert cdf_distance([1, 2], [2, 4]) == 0
    assert np.array_equal(rejection_indices([.1, .2, .20001], .2), [0, 1])


def test_read_stream_replay_and_fixed_depth():
    a = observe_counts([2, 0, 7], 1001, np.random.default_rng(42))
    b = observe_counts([2, 0, 7], 1001, np.random.default_rng(42))
    assert np.array_equal(a, b)
    assert a.sum() == 1001 and a[1] == 0


@pytest.mark.parametrize("bad", [[1.2, 3], [-1, 3], [np.nan, 1], [], [[1, 2]], [2**53, 0]])
def test_bad_population_counts_raise(bad):
    with pytest.raises(ValueError):
        observe_counts(bad, 10, np.random.default_rng(1))


def test_invalid_distance_inputs_raise():
    with pytest.raises(ValueError):
        cdf_distance([1, 2], [0, 0])
    with pytest.raises(ValueError):
        cdf_distance([1], [1, 2])
    for eps in [np.inf, np.nan, -1, 1.01]:
        with pytest.raises(ValueError):
            rejection_indices([.1, np.inf], eps)


def pilot_module():
    # The pilot is an executable example, deliberately outside the public API.
    import importlib.util
    from pathlib import Path
    import sys
    examples = Path(__file__).resolve().parents[1] / 'examples'
    sys.path.insert(0, str(examples))
    try:
        spec = importlib.util.spec_from_file_location('abc_pilot', examples / 'abc_pilot.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        return module
    finally:
        sys.path.remove(str(examples))


def test_pilot_bins_preserve_rounding_window_and_maximum(tmp_path):
    m = pilot_module()
    f = tmp_path / 'data.csv'
    f.write_text('id,chr,length,freq\na,chrM,39,7\na,chrM,40,3\na,chrM,1000,2\na,chrM,16569,1\n')
    labels, h = m.read_data(f, m.edges_bp())
    assert labels[0]['source_n'] == 13 and labels[0]['excluded'] == 7
    assert h.sum() == 6 and h[0,0] == 3 and h[0,-1] == 1
    assert np.cumsum(h[0])[np.where(m.edges_bp()[1:] == 1000.5)[0][0]] == 5


def test_pilot_resource_errors_propagate_without_retry():
    from types import SimpleNamespace
    from io import StringIO
    from frime_forest.trees import ResourceLimitError
    m = pilot_module()
    class Fail:
        calls = 0
        def run(self, *args, **kwargs):
            self.calls += 1
            raise ResourceLimitError('deliberate cap')
    runner = Fail()
    log = m.RunLog(runner, SimpleNamespace(seed=1,delta=1e-10), StringIO())
    with pytest.raises(ResourceLimitError, match='deliberate cap'):
        log.population(m.TRUTH, 300, 0)
    assert runner.calls == 1 and log.n == 0


def test_pilot_stage_streams_and_boundary_rate_gauge():
    from frime_forest import replicate_seed
    m = pilot_module()
    stages = list(m.STAGES.values())
    roots = [replicate_seed(20261012, s) for s in stages]
    assert len(set(roots)) == len(stages)
    a = m.prior_draws(10, 123, 300)
    b = m.prior_draws(10, 123, 400)
    assert not np.array_equal(a, b)
    model = m.model_for(m.TRUTH)
    boundary = 500 / m.L_BP
    assert model.rho(boundary) == pytest.approx(1)
    assert model.exit.rate(boundary/2, 1) == pytest.approx(1.5)
    assert model.c_i == 30


def test_marked_tree_projection_at_events_inherits_ancestor_pruning():
    from dataclasses import replace
    from frime_forest import Model, Exit
    from frime_forest.trees import make_skeleton, Skeleton
    base=Model(L=1,ell=.02,c_f=2,alpha=0,c_i=5,exit=Exit('constant',c_e=1))
    original=make_skeleton(base,0,723,root_size=1)
    cutoff=.15
    keep=np.flatnonzero(original.sizes>cutoff)
    remap={int(old):new for new,old in enumerate(keep)}
    parents=np.array([-1 if original.parents[k]<0 else remap[int(original.parents[k])] for k in keep])
    projected_model=replace(base,ell=cutoff)
    for exit_rate in (0.,1.,100.):
        model=replace(base,exit=Exit('constant',c_e=exit_rate))
        mproj=replace(projected_model,exit=model.exit)
        projected=Skeleton(original.family_id,original.root_size,mproj.fragmentation_key(),
            original.sizes[keep],original.births[keep],original.splits[keep],parents,original.exit_marks[keep])
        f=original.prune(model); g=projected.prune(mproj)
        times=np.unique(np.r_[0,f.births,f.ends])
        times=times[np.isfinite(times)]
        times=np.r_[times,(times[1:]+times[:-1])/2,times[-1]+1]
        for t in times:
            full=f.at(float(t))
            assert np.array_equal(full[full>cutoff],g.at(float(t)))
