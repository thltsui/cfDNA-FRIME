from dataclasses import replace
import math
import numpy as np
import pytest

from frime_forest import (Model, Exit, certify, Skeleton, make_family, make_skeleton,
                          make_plan, presample, assemble, sample, ResourceLimitError)
from frime_forest.certificate import contraction


@pytest.mark.parametrize("kwargs", [{"ell": 0}, {"ell": 1}, {"c_f": 0}, {"c_i": -1},
                                    {"a": 0}, {"alpha": math.nan}, {"L": math.inf}])
def test_invalid_model(kwargs):
    with pytest.raises((ValueError, FloatingPointError)):
        Model(**kwargs)


@pytest.mark.parametrize("kind", ["constant", "pnb", "pfb"])
def test_exit_laws(kind):
    exit = Exit(kind, c_e=2, power=-1, boundary=0.4)
    xs = np.array([0.1, 0.3, 0.4, 0.8, 1.0])
    expected = {"constant": np.full(5, 2.0), "pnb": 2/xs,
                "pfb": 2*np.maximum(1/xs - 1/0.4, 0)}[kind]
    np.testing.assert_allclose([exit.rate(float(x),1) for x in xs], expected)


@pytest.mark.parametrize("delta", [0, 1, -0.1, math.nan])
def test_invalid_delta(delta):
    with pytest.raises(ValueError):
        certify(Model(), delta)


def test_uniform_split_contraction():
    for p in (2,3,4,6,8):
        assert contraction(1,1,p) == pytest.approx(1-2/(p+1))


def test_certificate_formula_and_monotonicity():
    m = Model(alpha=0, ell=.1, c_f=3, c_i=10, exit=Exit("constant", c_e=2))
    loose, tight = certify(m, 1e-3, powers=(2,)), certify(m, 1e-9, powers=(2,))
    A = (1-.1**3)/(3*.1**2)
    assert math.exp(tight.log_A) == pytest.approx(A)
    assert tight.kappa == pytest.approx(3)
    assert tight.T == pytest.approx(math.log(10*A/(3*1e-9))/3)
    assert tight.T > loose.T
    assert tight.mismatch_bound <= 1e-9*(1+1e-10)
    assert certify(replace(m, c_i=100), 1e-9, powers=(2,)).T > tight.T


@pytest.mark.parametrize("alpha", [-2, 0, 1, 2])
@pytest.mark.parametrize("kind", ["constant", "pnb", "pfb"])
def test_interval_bound_is_lower_than_rate(alpha, kind):
    m = Model(alpha=alpha, exit=Exit(kind))
    cert = certify(m, powers=(2,))
    xs = np.geomspace(m.ell, m.L, 2000)
    values = [contraction(m.a,m.b,2)*m.rho(float(x))+m.exit.rate(float(x),m.L) for x in xs]
    assert cert.kappa <= min(values)*(1+1e-12)
    # Refining the interval lower envelope should not degrade the bound.
    assert cert.kappa >= certify(m, powers=(2,), intervals=1).kappa*(1-1e-12)


def manual_skeleton():
    m = Model(ell=.1, exit=Exit(c_e=1))
    tree = Skeleton(0, 1, m.fragmentation_key(), np.array([1,.4,.6]), np.array([0,2,2]),
                    np.array([2,5,5]), np.array([-1,0,0]), np.array([1,10,10]))
    return m, tree


def test_ancestor_exit_prunes_both_children():
    m, tree = manual_skeleton()
    f = tree.prune(m)
    np.testing.assert_array_equal(f.at(.5), [1])
    assert len(f.at(1)) == len(f.at(2)) == len(f.at(3)) == 0


def test_birth_inclusive_and_end_exclusive():
    m, tree = manual_skeleton()
    f = tree.prune(replace(m, exit=Exit(c_e=0)))
    np.testing.assert_array_equal(f.at(0), [1])
    np.testing.assert_array_equal(f.at(2), [.4,.6])
    assert len(f.at(5)) == 0


def test_prune_is_nonmutating_and_checks_compatibility():
    m, tree = manual_skeleton()
    before = tree.splits.copy()
    tree.prune(m)
    tree.prune(replace(m, exit=Exit("pnb")))
    np.testing.assert_array_equal(tree.splits, before)
    with pytest.raises(ValueError):
        tree.prune(replace(m, alpha=1))


def test_skeleton_sizes_mass_and_cutoff():
    m = Model(ell=.02)
    for i in range(50):
        tree = make_skeleton(m, i, 5, root_size=1)
        assert np.all(tree.sizes > m.ell)
        assert np.all(tree.parents < np.arange(len(tree.sizes)))
        for t in (0,.2,1,3,10):
            f = tree.prune(replace(m, exit=Exit(c_e=0)))
            assert f.at(t).sum() <= 1+1e-12
        for j in range(len(tree.sizes)):
            children = tree.sizes[tree.parents == j]
            assert children.sum() <= tree.sizes[j]+1e-12
            assert 0 <= len(children) <= 2


def test_root_at_cutoff_is_empty():
    m = Model()
    assert len(make_skeleton(m,0,5,root_size=m.ell).sizes) == 0


def test_no_immigration():
    f = sample(Model(c_i=0), horizon=5)
    assert f.plan.certificate.T == 0
    assert len(f.at(0)) == len(f.at(5)) == 0


def test_resource_caps_abort():
    m = Model(ell=.01, exit=Exit(c_e=0))
    with pytest.raises(ResourceLimitError):
        make_skeleton(m,0,1,root_size=1,max_nodes=1)
    with pytest.raises(ResourceLimitError):
        make_plan(m,max_families=1)
    with pytest.raises(ResourceLimitError):
        sample(m,max_total_nodes=1)


def test_bank_assembly_no_recycling():
    m = Model(ell=.2,c_i=1)
    plan = make_plan(m,horizon=1,seed=25)
    n = len(plan.arrivals)
    bank = presample(m,n,seed=25)
    f = assemble(plan,bank)
    assert len(f.families) == n
    np.testing.assert_array_equal(f.at(.5), f.at(.5))
    with pytest.raises(ValueError):
        assemble(plan,bank[:-1])
    if n > 1:
        with pytest.raises(ValueError):
            assemble(plan,(bank[0],)*n)


def test_query_limits_and_histogram():
    f = sample(Model(ell=.1,c_i=2),horizon=2,seed=4)
    for t in (-1,2.1,math.nan):
        with pytest.raises(ValueError):
            f.at(t)
    bins = np.linspace(.1,1,11)
    np.testing.assert_array_equal(f.histogram(bins,.5), np.histogram(f.at(.5),bins)[0])
    for bins in ([0,0,1], [0,math.nan], [1], [[0,1]]):
        with pytest.raises(ValueError):
            f.histogram(bins)


@pytest.mark.parametrize("mode", ["lazy", "skeleton"])
def test_worker_and_batch_invariance(mode):
    m = Model(ell=.1,c_i=2)
    serial = sample(m,horizon=.5,seed=2026,mode=mode,workers=1,batch_size=1)
    parallel = sample(m,horizon=.5,seed=2026,mode=mode,workers=2,batch_size=7)
    np.testing.assert_array_equal(serial.plan.arrivals,parallel.plan.arrivals)
    for t in (0,.1,.5):
        np.testing.assert_array_equal(serial.at(t),parallel.at(t))


def test_presampled_parallel_bank_matches_serial():
    m = Model(ell=.1)
    a = presample(m,12,seed=16,workers=1,batch_size=2)
    b = presample(m,12,seed=16,workers=2,batch_size=5)
    for s,t in zip(a,b):
        for name in ("sizes","births","splits","parents","exit_marks"):
            np.testing.assert_array_equal(getattr(s,name),getattr(t,name))


def test_realised_lifetimes_are_not_cut_off_at_T():
    f = sample(Model(ell=.1,c_i=2,exit=Exit(c_e=0)),horizon=2,seed=6,mode="skeleton")
    assert all(math.isinf(family.horizon) for family in f.families)
    # T is only present in the immigration plan; each family has its own clocks.
    assert all(s >= -f.plan.certificate.T for s in f.plan.arrivals)


def test_parallel_pruning_of_precomputed_bank():
    m = Model(ell=.1,c_i=2)
    plan = make_plan(m,seed=34)
    bank = presample(m,len(plan.arrivals),seed=34)
    one = assemble(plan,bank)
    two = assemble(plan,bank,workers=2,batch_size=3)
    np.testing.assert_array_equal(one.at(),two.at())


@pytest.mark.parametrize("workers", [0,-1,1.5,True])
def test_worker_validation(workers):
    with pytest.raises(ValueError):
        sample(Model(c_i=0),workers=workers)
