"""Seeded statistical checks against analytic identities, not legacy-code output.

Six-Monte-Carlo-standard-error thresholds are deliberately non-brittle. Passing
these tests is evidence, not a substitute for the construction proof.
"""
from dataclasses import replace
import math
import numpy as np
import pytest
from frime_forest import Model, Exit, make_skeleton, sample, certify


def assert_mc_mean(values, expected, sigma=6):
    values = np.asarray(values)
    se = values.std(axis=0,ddof=1)/math.sqrt(len(values))
    assert np.all(np.abs(values.mean(axis=0)-expected) <= sigma*se+1e-9)


def test_uniform_full_tree_expected_split_count():
    m = Model(ell=.1,exit=Exit(c_e=0))
    # Every stored pure-tree node eventually splits.
    counts = [len(make_skeleton(m,i,2026,root_size=1).sizes) for i in range(2000)]
    assert_mc_mean(counts, 2/m.ell-1)
    counts = [len(make_skeleton(m,i,101).sizes) for i in range(2000)]
    assert_mc_mean(counts, m.L/m.ell-1)


def test_root_competing_exponentials():
    m = Model(ell=.6,c_f=2,alpha=0,exit=Exit(c_e=3))
    holds, splits = [], []
    for i in range(5000):
        tree = make_skeleton(m,i,728,root_size=1)
        kill = tree.exit_marks[0]/3
        holds.append(min(tree.splits[0],kill))
        splits.append(tree.splits[0] < kill)
    assert_mc_mean(holds,1/5)
    assert_mc_mean(splits,2/5)


def test_uniform_constant_rate_stationary_histograms():
    # Derivation: h(x)=L^-1 (L/x)^(2c_f/(c_f+c_e)); mean density = c_i*h/(c_f+c_e).
    m = Model(ell=.2,c_f=2,alpha=0,c_i=5,exit=Exit(c_e=1))
    edges = np.linspace(m.ell,m.L,6)
    exponent = 2*m.c_f/(m.c_f+m.exit.c_e)
    scale = m.c_i/(m.L*(m.c_f+m.exit.c_e))*m.L**exponent
    expected = scale*(edges[1:]**(1-exponent)-edges[:-1]**(1-exponent))/(1-exponent)
    for mode,offset in (("lazy",0),("skeleton",10000)):
        h = [sample(m,seed=i+offset,mode=mode,delta=1e-9).histogram(edges) for i in range(400)]
        assert_mc_mean(h,expected)


def test_stationary_mean_without_exit():
    m = Model(ell=.25,c_f=2,alpha=0,c_i=3,exit=Exit(c_e=0))
    counts, masses = [], []
    for i in range(400):
        values = sample(m,seed=i+20321,delta=1e-9).at()
        counts.append(len(values)); masses.append(values.sum())
    assert_mc_mean(counts, m.c_i*m.L/m.c_f*(1/m.ell-1/m.L))
    assert_mc_mean(masses, m.c_i*m.L/m.c_f*math.log(m.L/m.ell))


def test_family_survival_and_p_mass_bound():
    m = Model(ell=.2,c_f=2,alpha=0,exit=Exit(c_e=.5))
    cert = certify(m,powers=(2,))
    ages = np.array([1,3,5])
    p_mass = []
    alive = []
    for i in range(1500):
        f = make_skeleton(m,i,345,root_size=1).prune(m)
        vals = [f.at(float(t)) for t in ages]
        p_mass.append([np.square(x).sum() for x in vals])
        alive.append([len(x)>0 for x in vals])
    p_mass = np.asarray(p_mass)
    bound = np.exp(-cert.kappa*ages)
    assert np.all(p_mass.mean(axis=0) <= bound+6*p_mass.std(axis=0)/math.sqrt(len(p_mass)))
    observed = np.mean(alive,axis=0)
    qbound = np.minimum(1,bound/m.ell**2)
    assert np.all(observed <= qbound+6*np.sqrt(np.maximum(observed*(1-observed),1/len(alive))/len(alive)))


def _chronological_family(model, root, horizon, rng):
    """Independent small CTMC reference, only for tests. No legacy imports."""
    sizes = [root] if root > model.ell else []
    t = 0.0
    while sizes:
        rho = np.array([model.rho(x) for x in sizes])
        exit_rates = np.array([model.exit.rate(x,model.L) for x in sizes])
        rates = np.concatenate((rho,exit_rates))
        total = rates.sum()
        t += rng.exponential(1/total)
        if t > horizon:
            break
        j = int(rng.choice(len(rates),p=rates/total))
        if j >= len(sizes):
            sizes.pop(j-len(sizes))
        else:
            x = sizes.pop(j)
            r = rng.beta(model.a,model.b)
            sizes.extend(y for y in (x*r,x*(1-r)) if y > model.ell)
    return np.array(sizes)


@pytest.mark.parametrize("alpha", [-1, 0, 1])
@pytest.mark.parametrize("kind", ["constant", "pnb", "pfb"])
def test_nonuniform_beta_against_independent_chronological_reference(alpha,kind):
    from frime_forest import make_family
    m = Model(ell=.1,alpha=alpha,a=2,b=5,exit=Exit(kind,power=-.5))
    rng = np.random.default_rng(1729)
    forest_stats, reference_stats = [], []
    def stats(x):
        return [len(x), x.sum(), len(x)**2]
    for i in range(600):
        f = make_family(m,i,7743,1.0)
        x = f.at(1.0)
        y = _chronological_family(m,float(rng.uniform(0,m.L)),1.0,rng)
        forest_stats.append(stats(x)); reference_stats.append(stats(y))
    a,b = np.asarray(forest_stats),np.asarray(reference_stats)
    se = np.sqrt(a.var(axis=0,ddof=1)/len(a)+b.var(axis=0,ddof=1)/len(b))
    assert np.all(np.abs(a.mean(axis=0)-b.mean(axis=0)) <= 6*se+1e-9)
