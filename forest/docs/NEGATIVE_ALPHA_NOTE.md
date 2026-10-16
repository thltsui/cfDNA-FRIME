# Negative fragmentation exponents, cutoffs, and stationary singularities

Research note — 14 October 2026. These calculations concern **uniform binary splitting**, total immigration intensity $C_I>0$ and $C_F>0$, and the **expected number density** $u(x)$; the expected mass density is $xu(x)$. Immigration and splitting Beta parameters are distinct. Hold all rates and the normalized immigrant law fixed when varying the cutoff; remapping $C_I$ to preserve an abundance prior would change the question.

For every fixed $\ell>0$, only $x>\ell$ is retained. Rates are bounded on that compact size interval. The current simulator/certificate concerns this model. A density singularity at zero is outside its domain and does not invalidate a positive-cutoff run. Removing the cutoff additionally requires a construction allowing infinitely many splitting events and a specified treatment of dust.

**Biological reconciliation: retain a positive destruction scale.** DNA length is discrete, so an unlimited cascade to arbitrarily small positive length is a continuum idealization. Treating fragments at or below a fixed $\ell$ as irreversibly destroyed gives a consistent retained-fragment model. For $\alpha<0$,

$$\sup_{\ell<x\le L}\rho(x)=C_F\ell^\alpha<\infty.$$

The PDE above $\ell$ retains its gain/loss form. Offspring at or below $\ell$ are absorbed; this is a jump-fragmentation loss, **not** a requirement that $u(\ell)=0$. For uniform splitting the mean mass balance is

$$\frac{d}{dt}\int_\ell^L xu_t(x)\,dx
=C_I\int_\ell^L xp_I(x)\,dx
-\int_\ell^L xE(x)u_t(x)\,dx-J_t(\ell),
\qquad J_t(\ell)=\ell^2\int_\ell^L\frac{\rho(y)u_t(y)}y\,dy.$$

At equilibrium, retained immigrant mass is balanced by exit and destruction. The appropriate inference problem keeps $\ell>0$ fixed; no $\ell\downarrow0$ limit is required. This removes the zero-size singularity obstruction but is not itself a proof of every hydrodynamic tightness assumption. The hard destruction threshold is a modeling assumption; discreteness alone does not identify a biological degradation threshold.

Current ABC documentation describes a 1 bp biological cutoff and computes the exact projection above 39.5 bp (rounded reads >=40 bp). Any lower destruction cutoff gives the same observable process under common immigration, rates and node marks, since descendants only shrink and do not interact. These data cannot identify a destruction scale below that window. The observable threshold must not be called a measured biological destruction length. A continuous-size model above cutoff is already sufficient for this reconciliation; implementing integer splits would change the split kernel and needs separate validation.

The thesis specifies $I(x)=C_Ip_I(x)$, with $p_I=1/L$ for uniform immigration: [FRIME_math.tex, lines 134–155](</Users/holungtsui/Documents/GitHub/DPhil-Thesis/OxThesis-master/text/FRIME_math.tex:134>). Its stationary balance is at [lines 454–476](</Users/holungtsui/Documents/GitHub/DPhil-Thesis/OxThesis-master/text/FRIME_math.tex:454>). Lines 491 and 509–519 contain inconsistent boundary/density normalizations; their prefactors should not be imported silently.

**New deductions from stationary balance.** Put $\rho(x)=C_Fx^\alpha$ and $h=\rho+E$. For $\ell<x<L$, formally extended to $0<x<L$ in the cutoff limit,

$$h(x)u(x)=I(x)+2\int_x^L\frac{\rho(y)u(y)}y\,dy.$$

Define the downward fragmentation mass flux

$$J(x)=x^2\int_x^L\frac{\rho(y)u(y)}y\,dy.$$

Direct differentiation gives

$$u=\frac{I+2J/x^2}{h},\qquad
J'=xEu-xI=\frac{2E}{xh}J-\frac{x\rho I}{h},\qquad J(L)=0.$$

This integral boundary condition also handles Beta densities singular at $L$; a finite pointwise value $u(L)$ need not exist. At positive cutoff, removed mass includes both $J(\ell)$ and immigrants born below $\ell$.

With **no exit**, $J(x)=C_I\int_x^L yp_I(y)\,dy$. Consequently uniform immigration on $(0,L)$ gives the exact formula

$$\boxed{u(x)=\frac{C_IL}{C_F}x^{-\alpha-2}},\qquad \ell<x<L.$$

For scaled Beta immigration with shapes $a_I,b_I>0$,

$$u(x)\sim\frac{2C_I\mathbb E[X]}{C_F}x^{-\alpha-2}\quad(x\downarrow0).$$

Thus $\int_0^L u<\infty$ exactly when $\alpha<-1$, whereas $\int_0^L xu<\infty$ exactly when $\alpha<0$. For $-2<\alpha<-1$ the density diverges but is integrable; at $\alpha=-2$ it approaches a constant; below $-2$ it vanishes. At $\alpha=-1$, expected number diverges logarithmically; for $-1<\alpha<0$, it diverges as a power while expected mass remains finite. These are **expectation/integrability conclusions**. Infinite expectation alone does not imply almost-sure infinity.

**New formal PFB asymptotics.** Take the current $\alpha_E=-1$ law

$$E(x)=k(1/x-1/B)_+,\qquad k=C_EL>0.$$

The thesis PFB discussion [lines 550–557](</Users/holungtsui/Documents/GitHub/DPhil-Thesis/OxThesis-master/text/FRIME_math.tex:550>) specializes to different fragmentation parameters; the following extensions are deductions, not quoted theorems. Write $p_I(x)\sim c x^{a_I-1}$; uniform immigration corresponds to $a_I=1$.

* If $\alpha<-1$, fragmentation dominates exit. Since $2E/(xh)$ is integrable at zero,
  $$J_0=\int_0^L\frac{y\rho(y)I(y)}{h(y)}
  \exp\!\left[-\int_0^y\frac{2E(z)}{zh(z)}\,dz\right]dy>0,$$
  and $u(x)\sim(2J_0/C_F)x^{-\alpha-2}$. The exponent is unchanged, but exit attenuates the coefficient.
* If $\alpha=-1$, let $\gamma=(k-C_F)/(k+C_F)\in(-1,1)$. For uniform immigration, $u\sim A x^\gamma$. For Beta immigration, the same holds when $a_I>\gamma$; equality gives $A x^{a_I}\log(L/x)$; $a_I<\gamma$ gives $A x^{a_I}$. The latter two cases require $\gamma>0$.
* If $-1<\alpha<0$, exit dominates: uniform immigration gives $u\sim Ax$. Beta immigration gives order $x^{\min(a_I,1)}$; when $a_I<1$, the coefficient is $C_Ic/k$.

Here $A>0$ depends on the full model. All these positive-$k$ PFB regimes have finite expected number and mass. The flux identity yields

$$J_0=C_I\mathbb E[X]-\int_0^L xE(x)u(x)\,dx.$$

Thus the formal zero-cutoff balance has positive dust flux for $\alpha<-1$ and zero dust flux for $-1\le\alpha<0$. This is a boundary-mass statement, not a proof of pathwise well-posedness for a new simulator.

**External theorem context.** [Haas (2005), §1.1 and Theorem 7](https://www.math.univ-paris13.fr/~haas/AAPFragmentationImmigration.pdf), *Annals of Applied Probability* 15, 1958–1996, treats self-similar fragmentation with immigration. For bounded immigrant sizes, its negative-index results give stationary states with finite mass almost surely; for $\alpha<-1$ they have finitely many positive fragments almost surely. Section 1.1 discusses finite-time extinction to dust. These results support the no-exit interpretation; they are not citations for the new PFB asymptotics above.

Keep three limits separate: $n\to\infty$ averages independent populations; $t\to\infty$ concerns stationarity; $\ell\downarrow0$ changes the retained state space. The thesis [lines 245–292 and 356–369](</Users/holungtsui/Documents/GitHub/DPhil-Thesis/OxThesis-master/text/FRIME_math.tex:245>) develops the mean/hydrodynamic limit, with tightness assumed. Count/mass laws of large numbers require their respective finite first moments. Exchanging limits requires justification.

The current history bound is not uniform in $\ell$: for negative $\alpha$, $\kappa_p=(1-m_p)C_FL^\alpha$ for PFB, while $A_p\asymp\ell^{-p}$. Its sufficient lookback therefore grows like $(p/\kappa_p)\log(L/\ell)$. This does not prove actual extinction times diverge. Separately, monotone fragmentation makes restriction above any fixed observable threshold exact: discarded small fragments cannot return. When a zero-cutoff coupling has integrable $u$, its spatial-cutoff mismatch probability is at most $\int_0^\ell u(x)dx$; finite mass alone only controls omitted mass.

Current S04 candidate 3269 has $\alpha=-0.9602$, while 3243 and 3034 have $\alpha=-1.3445,-1.4803$. Their fresh median CDF gaps are approximately $0.0303,0.0453,0.0451$ ([saved results](</Users/holungtsui/Documents/GitHub/cfDNA-FRIME/forest/results/immigration-sensitivity-20261014/full-screen/confirmation-summary.json>)). Acceptable observed-window profiles can therefore occur on both sides of $-1$. These comparisons do not identify the zero-size boundary mechanism, an absolute biological timescale, or a unique fragmentation exponent.
