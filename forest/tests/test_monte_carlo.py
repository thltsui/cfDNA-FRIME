"""Population-level MC checks: law, independent IDs, summary reductions, failures."""
from dataclasses import replace
import math

import numpy as np
import pytest

from frime_forest import (Model, Exit, MonteCarlo, monte_carlo, replicate_seed,
                          sample, ResourceLimitError)


@pytest.mark.parametrize('mode', ['lazy', 'skeleton'])
def test_mc_matches_existing_full_population_sampler(mode):
    m = Model(L=2, ell=.2, alpha=.5, a=2, b=5, c_i=3, exit=Exit('pfb', boundary=.8))
    bins = np.linspace(.3, 1.7, 9)  # A subwindow must not redefine total count/mass.
    out = monte_carlo(m, 5, bins=bins, seed=73, start_replicate=8, mode=mode)
    for row, seed in enumerate(out.seeds):
        f = sample(m, seed=seed, mode=mode)
        values = f.at()
        np.testing.assert_array_equal(out.histograms[row], np.histogram(values, bins)[0])
        assert out.counts[row] == len(values)
        assert out.masses[row] == math.fsum(float(family.at(-float(s)).sum())
                                           for s, family in zip(f.plan.arrivals, f.families))
        assert out.generated_nodes[row] == f.stored_nodes
        assert out.immigrants[row] == len(f.families)
        assert out.counts[row] >= out.histograms[row].sum()


@pytest.mark.parametrize('mode', ['lazy', 'skeleton'])
def test_mc_ids_chunking_and_persistent_worker_batch_invariance(mode):
    m = Model(ell=.1, c_i=3, alpha=-1, a=.5, b=.5, exit=Exit('pnb'))
    bins = np.linspace(.1, 1, 7)
    serial = monte_carlo(m, 6, bins=bins, seed=998, mode=mode, batch_size=1)
    with MonteCarlo(workers=2, batch_size=7) as runner:
        pool = runner._pool
        first = runner.run(m, 2, bins=bins, seed=998, mode=mode)
        # An unrelated previous call cannot affect the streams of a later call.
        runner.run(replace(m, exit=Exit(c_e=0)), 1, bins=bins, seed=12)
        rest = runner.run(m, 4, bins=bins, seed=998, start_replicate=2, mode=mode)
        assert runner._pool is pool
        reversed_rows = [runner.run(m, 1, bins=bins, seed=998, start_replicate=i, mode=mode)
                         for i in range(5, -1, -1)]
    for field in ('replicate_ids', 'histograms', 'counts', 'masses', 'immigrants', 'generated_nodes'):
        combined = np.concatenate([getattr(first, field), getattr(rest, field)])
        np.testing.assert_array_equal(getattr(serial, field), combined)
        reordered = np.concatenate([getattr(r, field) for r in reversed_rows])[::-1]
        np.testing.assert_array_equal(getattr(serial, field), reordered)
    assert serial.seeds == first.seeds + rest.seeds
    assert len(set(serial.seeds)) == 6
    assert runner._pool is None
    with pytest.raises(RuntimeError, match='with block'):
        runner.run(m, 1, bins=bins)


def test_replicate_seed_encoding_does_not_overlap_adjacent_ensembles():
    seeds = {replicate_seed(root, i) for root in range(20) for i in range(20)}
    assert len(seeds) == 400
    assert replicate_seed(2**140, 2**50) == replicate_seed(2**140, 2**50)
    for value in (-1, True, 1.5):
        with pytest.raises(ValueError):
            replicate_seed(value, 0)
        with pytest.raises(ValueError):
            replicate_seed(0, value)


def test_mc_empty_population_mcse_and_joint_history_tolerance():
    out = monte_carlo(Model(c_i=0), 3, bins=[.1, .5, 1])
    assert not out.counts.any() and not out.histograms.any() and not out.masses.any()
    assert out.joint_history_bound == 0
    assert out.statistics()['count'] == {'mean': 0., 'mcse': 0.}
    one = monte_carlo(Model(c_i=1), 1, bins=[0,1], delta=1e-5)
    assert one.statistics()['histogram']['mcse'] is None
    many = monte_carlo(Model(c_i=1), 4, bins=[0,1], delta=1e-5)
    assert one.certificate == many.certificate  # delta is per replicate, not delta/n.
    assert many.joint_history_bound == pytest.approx(-math.expm1(-4*many.certificate.lambda_bound))
    assert many.joint_history_bound <= 4e-5
    with pytest.raises(ValueError):
        many.histograms[0,0] = 0


@pytest.mark.parametrize('workers', [1,2])
@pytest.mark.parametrize('cap', ['max_families', 'max_nodes', 'max_total_nodes'])
def test_mc_cap_failure_aborts_with_replicate_identity(workers, cap):
    m = Model(ell=.01, c_i=5, exit=Exit(c_e=0))
    with MonteCarlo(workers=workers, **{cap:1}) as runner:
        with pytest.raises(ResourceLimitError, match='replicate 17:'):
            runner.run(m, 3, bins=[0,1], seed=77, start_replicate=17)
        with pytest.raises(RuntimeError, match='with block'):
            runner.run(m, 1, bins=[0,1])


@pytest.mark.parametrize('kwargs', [{'n':0}, {'n':True}, {'start_replicate':-1},
                                  {'start_replicate':2**63}, {'bins':[1,1]},
                                  {'bins':[0,math.inf]}, {'mode':'recycle'}])
def test_mc_input_validation(kwargs):
    with pytest.raises(ValueError):
        monte_carlo(Model(), **{'n':1, 'bins':[0,1], **kwargs})


def test_pfb_positive_exit_underflow_raises_without_changing_exact_zero():
    e = Exit('pfb', c_e=5e-324, power=-1, boundary=.4)
    with pytest.raises(FloatingPointError):
        e.rate(.399, 1)
    assert e.rate(.4, 1) == 0 and e.rate(.9,1) == 0
    assert Exit('pfb', c_e=0).rate(.399,1) == 0


@pytest.mark.parametrize('mode', ['lazy','skeleton'])
def test_mc_stationary_first_moments_with_total_immigration(mode):
    # L != 1 detects confusion between total C_I and immigration density C_I.
    m = Model(L=2, ell=.4, c_f=2, alpha=0, c_i=5, exit=Exit(c_e=1))
    bins = np.linspace(m.ell, m.L, 6)
    out = monte_carlo(m, 600, bins=bins, seed=691, mode=mode, delta=1e-9)
    exponent = 2*m.c_f/(m.c_f+m.exit.c_e)
    scale = m.c_i/(m.L*(m.c_f+m.exit.c_e))*m.L**exponent
    expected = scale*np.diff(bins**(1-exponent))/(1-exponent)
    summary = out.statistics()
    means = np.array(summary['histogram']['mean'])
    se = np.array(summary['histogram']['mcse'])
    assert np.all(np.abs(means-expected) <= 6*se+1e-9)
    np.testing.assert_allclose(se, out.histograms.std(axis=0,ddof=1)/math.sqrt(600))
    expected_mass = scale*(m.L**(2-exponent)-m.ell**(2-exponent))/(2-exponent)
    assert abs(summary['mass']['mean']-expected_mass) <= 6*summary['mass']['mcse']+1e-9


def test_mc_statistics_overflow_raises():
    out = monte_carlo(Model(c_i=0), 2, bins=[0,1])
    extreme = replace(out, masses=np.array([1e200, 2e200]))
    with pytest.raises(FloatingPointError):
        extreme.statistics()


def test_mc_context_misuse_raises():
    runner = MonteCarlo()
    with pytest.raises(RuntimeError):
        runner.run(Model(), 1, bins=[0,1])
    with runner:
        with pytest.raises(RuntimeError, match='already active'):
            runner.__enter__()
