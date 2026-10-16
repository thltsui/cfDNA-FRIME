"""Deterministic observable-abundance prior map, not a population sampler.

For uniform splitting, K=E_X integral_0^infinity |C_X(t)|dt and E|U|=c_i*K.
Numerical ODE/quadrature error affects this prior coordinate only. The history
certificate and all random populations continue to use the separate exact-law
family engine. SciPy is an optional study dependency, imported only on use.
"""
from dataclasses import dataclass, replace
import math

from .model import Model, finite


@dataclass(frozen=True)
class AbundanceCoefficient:
    value: float
    quadrature_error: float
    refinement_relative_difference: float
    requested_rtol: float
    ode_evaluations: int
    scipy_version: str
    final_rtol: float
    refinement_history: tuple[float, ...]
    method: str = 'DOP853 in log length; endpoint-weighted QUADPACK; tenfold refinement'


def _integrate(model, rtol):
    from scipy.integrate import solve_ivp, quad, IntegrationWarning
    from scipy.special import betaln
    import warnings
    d=model.ell/model.L
    scale=model.rho(model.L)
    span=math.log(model.L/model.ell)
    boundaries=[d,1.]
    if model.exit.kind=='pfb' and model.ell<model.exit.boundary<model.L:
        boundaries.insert(1,model.exit.boundary/model.L)

    def rates(z):
        rho=math.exp(model.alpha*math.log(z))
        raw_exit=model.exit.rate(model.L*z,model.L)
        exit=raw_exit/scale
        total=rho+exit
        if not math.isfinite(total) or rho<=0 or (raw_exit>0 and exit<=0):
            raise FloatingPointError('scaled occupancy rates are not representable')
        return rho,exit,total

    def rhs(t,w):
        z=math.exp(t-span)
        rho,exit,total=rates(z)
        return [1/(z*total)-2*(exit/total)*w[0]]

    pieces=[];initial=[0.];evaluations=0
    for lo,hi in zip(boundaries[:-1],boundaries[1:]):
        sol=solve_ivp(rhs,(math.log(lo/d),math.log(hi/d)),initial,
                      method='DOP853',rtol=rtol,atol=rtol*1e-3,dense_output=True)
        if not sol.success or not math.isfinite(float(sol.y[0,-1])):
            raise FloatingPointError('occupancy ODE failed: '+sol.message)
        pieces.append(sol);initial=sol.y[:,-1];evaluations+=sol.nfev
    ai,bi=(1.,1.) if model.immigration.is_uniform else (model.immigration.a,model.immigration.b)
    log_beta=float(betaln(ai,bi))
    values=[];errors=[]
    with warnings.catch_warnings():
        warnings.simplefilter('error',IntegrationWarning)
        for i,(lo,hi,sol) in enumerate(zip(boundaries[:-1],boundaries[1:],pieces)):
            def smooth(z):
                # Endpoint right limits are valid measure-zero extensions of v.
                w=float(sol.sol(math.log(z/d))[0])
                rho,exit,total=rates(z)
                v=(1+2*rho*z*w)/total
                answer=math.exp((ai-1)*math.log(z)-log_beta)*v
                if not math.isfinite(answer) or answer<=0:
                    raise FloatingPointError('positive occupancy integrand is not representable')
                return answer
            if hi==1.:
                value,error=quad(smooth,lo,hi,weight='alg',wvar=(0.,bi-1),
                                 epsabs=rtol*1e-3,epsrel=rtol,limit=200)
            else:
                value,error=quad(lambda z:smooth(z)*math.exp((bi-1)*math.log1p(-z)),
                                 lo,hi,epsabs=rtol*1e-3,epsrel=rtol,limit=200)
            values.append(value);errors.append(error)
    value=math.fsum(values)/scale
    error=math.fsum(errors)/scale
    if not math.isfinite(value) or value<=0 or not math.isfinite(error):
        raise FloatingPointError('positive occupancy coefficient is not representable')
    return value,error,evaluations


def abundance_coefficient(model: Model, *, rtol=2e-8) -> AbundanceCoefficient:
    """K for the model's retained window; c_i is intentionally ignored.

    The solver is checked by a tenfold tolerance refinement for every call.
    These numerical error diagnostics are not rigorous interval bounds.
    """
    import scipy
    if model.a!=1 or model.b!=1:
        raise ValueError('abundance prior mapping currently requires uniform splitting')
    finite(rtol,'rtol',positive=True)
    if not 1e-12<=rtol<=1e-4:
        raise ValueError('rtol must be between 1e-12 and 1e-4')
    coarse,_,evaluations=_integrate(model,rtol)
    history=[]
    for step in range(1,5):
        fine_rtol=rtol/(10**step)
        if fine_rtol<1e-13:
            break
        fine,error,n=_integrate(model,fine_rtol)
        evaluations+=n
        discrepancy=abs(coarse-fine)/fine
        history.append(discrepancy)
        if discrepancy<=10*rtol and error/fine<=10*rtol:
            return AbundanceCoefficient(fine,error,discrepancy,rtol,evaluations,
                                        scipy.__version__,fine_rtol,tuple(history))
        coarse=fine
    raise FloatingPointError('abundance mapping failed numerical refinement: '+repr(history))


def with_expected_abundance(model: Model, expected_count: float, *, rtol=2e-8):
    """Return a model with c_i=expected_count/K and recorded mapping diagnostics."""
    target=finite(expected_count,'expected_count',positive=True)
    coefficient=abundance_coefficient(model,rtol=rtol)
    intensity=target/coefficient.value
    if not math.isfinite(intensity) or intensity<=0:
        raise FloatingPointError('mapped immigration intensity is not representable')
    return replace(model,c_i=intensity),coefficient
