# FRIME as a Poisson forest: mathematical specification

**Status:** implementation specification and derived proofs for review, 8 September 2026.
This is not a claim that every formula below already appears in the thesis.

## 1. Source and scope

The model is Terence Tsui Ho Lung, *Branching and Fragmenting Particle Systems in
Biology*, Oxford DPhil thesis (17 November 2023): Definitions 7.1–7.4 (printed
pp. 143–145), Lemma 7.5 / equation (7.6) (p. 146), and Section 8.1.2 / equation
(8.2) (pp. 160–161). The thesis supplies binary Beta fragmentation, inherited
red/blue exit colouring, independent immigrant-family superposition, and the
proposal to truncate past immigration using a survival-tail estimate.

Sections 3–7 below spell out the implementation and derive an explicit moment
bound, a path-coupling guarantee, and analytic validation identities. These are
extensions worked out for this implementation, not quotations from the thesis.

Use the **Chapter 7 convention**: `c_i` is TOTAL immigration intensity. Equation
(7.5) gives intensity `c_i p_I(x) dx dt`, with `p_I=1/L` here. Chapter 6's algorithm
elsewhere uses `c_i L`; do not transfer numerical rate settings between conventions
without conversion. No old programme behaviour overrides this specification.

## 2. Implemented model and assumptions

For `0 < ell < L < infinity`, retain fragments of size `x > ell`. A fragment splits
at rate `rho(x)=c_f x^alpha`, `c_f>0`, into `xr` and `x(1-r)` with independent
`r~Beta(a,b)`, `a,b>0`. Descendants evolve independently. Immigration is a Poisson
process with rate `c_i>=0`, independently marked by sizes uniform on `[0,L]`.
Sizes never grow and there is no interaction between fragments or families.

Exit rate `E(x)` is one of:

- CON: `c_e`, allowing `c_e=0` as the no-exit control.
- PNB: `c_e (x/L)^beta`, `beta<0`.
- PFB: `c_e max((x/L)^beta-(B/L)^beta,0)`, `beta<0`, `0<B<=L`.

The thesis normalisation for PNB/PFB is `c_e=1`; an explicit nonnegative amplitude
is a documented extension. Instant removal at `ell` remains part of the model,
including when the ordinary exit rate is zero. Exponential clocks are in one
consistent time unit; `c_f` has units length^(-alpha)/time.

Rate functions are bounded on the retained compact size interval. A family's live
fragment count cannot exceed its root mass divided by ell. Its total event rate
is consequently bounded, so there is no finite-time explosion above the cutoff.

## 3. Marked-tree construction and equality in law

Each potential node v carries a length x_v, its parent, birth age b_v, a split
ratio R_v, and independent Exp(1) marks Z_v and Y_v. Its potential splitting time
is `s_v=b_v+Z_v/rho(x_v)` and its killing time is
`k_v=b_v+Y_v/E(x_v)`, interpreted as infinity when E=0.

If `k_v < s_v`, v exits and **every potential descendant is unreachable**.
Otherwise it splits and its children are born at s_v; retain only children above
ell. Shared ancestral killing is essential. Independent thinning of terminal
fragments has a different joint law.

Inductively, each reachable node has independent fragmentation and exit clocks,
with precisely the rates in Definition 7.2, and conditional on splitting its
children have fresh independent marks. Hence the reachable tree has the FRIME
family law. A local competing-clock implementation is equivalent: hold an
Exp(rho+E) duration, then split with probability rho/(rho+E). Do not combine that
Bernoulli probability with an independent Exp(rho) holding duration.

The pure skeleton is generated in depth-first, parent-before-child order.
Generating all trees before immigration times does not alter their joint law,
provided trees and immigration marks are independent. Each bank tree must be used
once within a population. Reusing a bank for paired parameter comparisons is
allowed but creates correlated outputs, not independent replicates.

A fragment is present at age u precisely when its ancestors split before exit
and `b_v <= u < min(s_v,k_v)`. `Forest.at(t)` is deterministic and never draws new
randomness. Both full and lazy traversals implement the same law; their seeded
realisations can differ because pruning changes the sequence of visited nodes.

## 4. Stationary construction and history error

Let C_x(u) be the surviving-family counting measure and tau_x its extinction age.
With Poisson immigration points `(S_i,X_i)` over the real line, write

$$U(t)=\sum_{S_i\leq t} C_{X_i}^{(i)}(t-S_i).$$

For a snapshot at 0 and lookback T, simulate `N~Poisson(c_i T)` independent trees,
with independent ages uniform on `[0,T]`. Their superposition is U_T. For a path
on `[0,H]`, instead use immigration times uniform on `[-T,H]` with count
`Poisson(c_i(T+H))` and natural family lifetimes.

By independent marking/thinning of the Poisson process, the number of omitted
families still alive at 0 is Poisson with mean

$$\Lambda_T=c_i\int_T^\infty \mathbb E_X[\Pr(\tau_X>u)]\,du
             =c_i\,\mathbb E[(\tau_X-T)_+].$$

Under the common-driver coupling,

$$\Pr(U_T(0)\ne U(0))=1-e^{-\Lambda_T},\qquad
 d_{TV}(\mathcal L(U_T(0)),\mathcal L(U(0)))\leq 1-e^{-\Lambda_T}.$$

On the event that no omitted family survives at 0, the two processes agree for
all future times if driven by the same later immigration. There is no resurrection.
Thus the same bound controls the path law on `[0,H]`, with no additional factor H.
**T truncates history, not family age.** Forcing every family to disappear at age
T would define a different finite-memory process and would need a different
path-error bound. This implementation does not do that.

## 5. Explicit non-asymptotic bound (derived here)

For p>1, set

$$V_p(C)=\sum_{x\in C}x^p,\qquad
m_p=\mathbb E[R^p+(1-R)^p]
=\frac{B(a+p,b)+B(a,b+p)}{B(a,b)}<1.$$

Dropping small children only increases p-mass loss. The family generator obeys

$$\mathcal G V_p(C)\leq-\sum_{x\in C}[(1-m_p)\rho(x)+E(x)]x^p.$$

Let kappa_p be any strictly positive lower bound for `(1-m_p)rho(x)+E(x)` over
`ell<x<=L`. The bounded-rate, bounded-state-moment setting justifies applying
Dynkin's formula and Gronwall to obtain

$$\mathbb E V_p(C_x(u))\leq x^p e^{-\kappa_pu},\qquad
\Pr(\tau_x>u)\leq \min\{1,(x/\ell)^p e^{-\kappa_pu}\}.$$

For x<=ell the family is identically empty. Integrating the tail proves finite
expected extinction time. Together with non-explosion, this also guarantees that
a size-truncated pure tree has finitely many nodes almost surely (apply the bound
with E=0). This is not a deterministic upper bound on realised node count.

For uniform X, define

$$A_p=\mathbb E[(X/\ell)^p 1_{X>\ell}]
=\frac{L^{p+1}-\ell^{p+1}}{(p+1)L\ell^p}.$$

Then

$$\Lambda_T\leq\frac{c_i A_p}{\kappa_p}e^{-\kappa_p T},\quad
T_\delta=\max\left\{0,\frac{1}{\kappa_p}
\log\frac{c_i A_p}{\kappa_p\delta}\right\}$$

gives `Lambda_T<=delta`, hence a mismatch bound at most delta. For c_i=0 use T=0.
A tenfold error reduction adds `log(10)/kappa_p` to this sufficient lookback.

### How the programme obtains a valid lower bound

Partition `[ell,L]` into deterministic intervals `[l_j,u_j]`. Since rho is
monotone and the supported exit functions are nonincreasing, use

$$k_j=(1-m_p)\min\{\rho(l_j),\rho(u_j)\}+E(u_j),\qquad
\kappa_p=\min_j k_j.$$

Every x in its interval satisfies the required inequality. **This is an interval
lower envelope, not a grid-point estimate of the minimum of the sum.** Finer
nested partitions can sharpen it without changing the proof. The implementation
checks integer moments p=2,3,4,6,8 and uses the shortest resulting bound; it does
not claim a globally optimal choice of p or T.

For p=2, `1-m_2=2ab/((a+b)(a+b+1))`; for uniform splits `m_p=2/(p+1)`.
Logs avoid large powers in A_p. A small roundoff margin is applied to kappa.

### Meaning and limits of a certificate

This is a bound on **past-history truncation for the ideal model**, evaluated
numerically in float64. It is not formal interval arithmetic and not a guarantee
on finite-precision/PRNG error, model misspecification, the choice of ell, or
Monte Carlo error in estimated moments. A theoretical TV guarantee against a
continuous law cannot literally apply to a discretely represented PRNG output.
Extremely small tolerances and extreme parameters may be unrepresentable; the
programme raises rather than silently changing rates or truncating trees.
Resource-limit failures are aborted attempts, not valid partial samples. Do not
repeatedly reject over-budget realisations and treat successful runs as iid.

The Haas tail expression used in thesis Section 8.1.2 is not used as an unqualified
all-time numerical bound in this implementation; no unknown asymptotic constants
are treated as 1.

## 6. Analytic validation identities (derived here)

For uniform splitting, the expected number n(x) of split nodes in the full pure
tree satisfies `n(x)=0` for x<=ell and

$$n(x)=1+\frac2x\int_\ell^x n(y)dy=2x/\ell-1\quad(x>\ell).$$

Thus for uniform immigration `E n(X)=L/ell-1`. Full-bank work depends strongly
on L/ell, even if pruning would eliminate most descendants.

For uniform splitting, alpha=0, and constant exit c_e, let
`s=c_f/(c_f+c_e)`. The expected density h of nodes actually born in one family,
averaged over uniform roots, solves
`h(x)=1/L+2s integral_x^L h(y)/y dy`. Therefore

$$h(x)=\frac1L(L/x)^{2s},\qquad
\mu(x)=\frac{c_i}{L(c_f+c_e)}(L/x)^{2s},\quad\ell<x<L.$$

Integrating mu over bins gives an exact stationary first-moment reference.
For c_e=0 the expected count is `c_i L/c_f (1/ell-1/L)` and expected retained mass
is `c_i L/c_f log(L/ell)`. The tests compare sample moments against these formulas.
A separate chronological CTMC reference checks nonuniform Beta splits, all three
exit families, and positive/negative alpha, including the second moment of count.
Six-standard-error statistical tests provide evidence, not proofs of identity.

## 7. Explicitly deferred extensions

A clock-integrated family-occupation estimator can estimate the stationary mean
using weights `w_v=product_ancestors rho/(rho+E)` and contributions
`w_v*g(x_v)/(rho(x_v)+E(x_v))`. It is **not** a stationary random realisation.
It is documented as a future task, not exposed as an implemented API here.

Dimensionless size skeletons can also be retimed: x_v=x*q_v gives
`D_v=x^(-alpha)/c_f * Z_v/q_v^alpha`. This scaling is not yet an exposed retiming
operation; v0.1 banks require the same fragmentation parameters. Exit/c_i may
change. Perfect infinite-past sampling, GPU/compiled kernels, clinical fitting,
unbounded immigrant sizes and interacting exit mechanisms are out of scope.
