from dataclasses import replace
import math
import numpy as np
import pytest

pytest.importorskip('scipy')
from scipy.special import betaincc
from frime_forest import Model,Exit,Immigration
from frime_forest.abundance_prior import abundance_coefficient,with_expected_abundance


@pytest.mark.parametrize('alpha',[-1.5,-1,-.3,0,1,1.5])
def test_uniform_no_exit_exact_integrated_family_count(alpha):
    m=Model(L=3,ell=.13,c_f=2.1,alpha=alpha,exit=Exit(c_e=0))
    power=alpha+1
    exact=3**(-alpha)/2.1*(math.expm1(power*math.log(3/.13))/power if power else math.log(3/.13))
    result=abundance_coefficient(m)
    assert result.value==pytest.approx(exact,rel=2e-8)


@pytest.mark.parametrize('a,b',[(.5,.5),(.5,8),(8,.5),(2,1),(1,2),(1,1)])
def test_beta_zero_power_no_exit_closed_form(a,b):
    m=Model(L=2,ell=.17,c_f=1.7,alpha=0,exit=Exit(c_e=0),immigration=Immigration('beta',a,b))
    q=m.ell/m.L
    exact=((2/q)*(a/(a+b))*betaincc(a+1,b,q)-betaincc(a,b,q))/m.c_f
    assert abundance_coefficient(m).value==pytest.approx(exact,rel=2e-8)


@pytest.mark.parametrize('a,b,sign',[(1,1,0),(2,1,1),(1,2,-1)])
def test_constant_exit_closed_form(a,b,sign):
    m=Model(L=3,ell=.2,c_f=1,alpha=0,exit=Exit(c_e=1),immigration=Immigration('beta',a,b))
    exact=.5*math.log(3/.2)+sign*(1-(.2/3)**2)/4
    assert abundance_coefficient(m).value==pytest.approx(exact,rel=2e-8)


def test_pfb_scaling_and_uniform_canonical_mapping():
    m=Model(L=2,ell=.02,c_f=2,alpha=-1.2,exit=Exit('pfb',c_e=.4,boundary=.3),
            immigration=Immigration('beta',.5,8))
    k=abundance_coefficient(m)
    scaled=replace(m,c_f=m.c_f*7,exit=replace(m.exit,c_e=m.exit.c_e*7))
    assert abundance_coefficient(scaled).value==pytest.approx(k.value/7,rel=2e-8)
    mapped,c=with_expected_abundance(m,1234)
    assert mapped.c_i*c.value==pytest.approx(1234)
    assert replace(mapped,c_i=m.c_i)==m
    u=replace(m,immigration=Immigration())
    assert abundance_coefficient(u)==abundance_coefficient(replace(u,immigration=Immigration('beta')))


def test_nonuniform_splitting_and_invalid_target_refused():
    with pytest.raises(ValueError,match='uniform splitting'):
        abundance_coefficient(Model(a=2))
    for x in (0,-1,np.inf,np.nan):
        with pytest.raises(ValueError):
            with_expected_abundance(Model(),x)
