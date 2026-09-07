# Wave 34 route C: resonant binary cost, moving-mass correction, and posterior exclusion

## Status

- **Verified:** every coordinate contribution to the weighted screened energy
  has an exact binary exponential-family representation.  With the edge
  marginal artificially held fixed, its integral is exactly a forward
  Bernoulli KL divergence, with all orientation and factor-of-two conventions
  retained.
- **Verified:** the actual edge marginal moves with `t`.  The difference
  between the weighted energy and the endpoint conditional KL is exactly one
  explicit covariance integral.  This is the missing term in any attempt to
  replace the resonant integral by endpoint binary divergences.
- **Verified and signing-specific:** on a vertex edge, the response is exactly
  the log ratio of posterior probabilities that the selector omits that
  vertex.  A resonant sign reversal therefore forces exponentially strong
  posterior inclusion at the endpoint which reverses the base Gibbs bias.
  This identity has no analogue for the global orientation edge, which must
  be kept separately.
- **Verified:** the best normalized scalar restoring fit has a closed form in
  terms of the angle between `X_t` and `B_t`.  This sharpens the diagnostic
  version of (10.961), but supplies no minimizer-specific lower bound by
  itself.
- **Open target:** a uniform exact-minimizer theorem controlling the endpoint
  conditional binary cost, its negative migration covariance, and the
  restoring angle/rate would prove the harmonic parent target with exactly
  the weight in (10.900).
- **Numerical finite walls:** exhaustive calculations on the exact order-nine
  minimizer show that resonant edges need not be a small part of the energy,
  the migration covariance need not have the favorable sign, and the small
  `A_4/A_6` restoring constants from Wave 33 are not stable across known exact
  minimizers.  These are finitary mechanism walls, not scalable
  counterexamples.

The exhaustive checker is `tmp/harmonic_resonance_r34_check.py`.  It uses the
repository virtual environment and enumerates every oriented cut, selector,
and coordinate edge.  Reported integrals and decimals are numerical
Gauss--Legendre evaluations; the displayed algebraic identities are proved
below independently of those evaluations.

## 1. Exact fixed-edge binary KL integral

Use the oriented-cut cube with coordinate involutions
`T_0,...,T_(n-1)`, where `T_0` is global orientation reversal and the
remaining coordinates are projective vertex flips.  Put

\[
 \mu_t(d)=\frac{\nu_\beta(d)e^{t g(d)}}{Z_t},\qquad
 g=\log f_\beta,\qquad 0\le t\le1.
\]

Fix one coordinate and one unordered edge `e={x,y}`, oriented temporarily
from `x` to `y`.  Define

\[
 \omega_e=\log\frac{\nu_\beta(y)}{\nu_\beta(x)},\quad
 \chi_e=g(y)-g(x),\quad
 u_e(t)=\omega_e+t\chi_e,
\]

\[
 M_e(t)=\mu_t(x)+\mu_t(y),\qquad
 p_e(t)=\frac{\mu_t(y)}{M_e(t)}
       =\frac{e^{u_e(t)}}{1+e^{u_e(t)}}.
\]

The conditional log partition, normalized at `t=0`, is

\[
 A_e(t)=\log\frac{\nu(x)e^{t g(x)}+\nu(y)e^{t g(y)}}
                    {\nu(x)+\nu(y)}.
\]

Direct differentiation gives

\[
 A'_e(t)=(1-p_e(t))g(x)+p_e(t)g(y),\qquad
 A''_e(t)=p_e(t)(1-p_e(t))\chi_e^2
 =\frac{\chi_e^2}{4\cosh^2(u_e(t)/2)}.                 \tag{R34.1}
\]

The unordered heat-bath conductance is

\[
 c_e(t)=M_e(t)p_e(t)(1-p_e(t)),
\]

so `c_e(t) chi_e^2=M_e(t)A''_e(t)`, exactly recovering
(10.957).  Now set

\[
 k_e(t)=tA'_e(t)-A_e(t).
\]

The conditional likelihood with respect to the base Bernoulli law is
`exp(tg-A_e(t))`; hence

\[
 \boxed{
 k_e(t)=D(\operatorname{Ber}(p_e(t))\Vert
          \operatorname{Ber}(p_e(0))),\qquad
 k'_e(t)=tA''_e(t).
 }                                                        \tag{R34.2}
\]

In particular the exact fixed-edge weighted integral is

\[
 \boxed{
 K_e:=\int_0^1 t\frac{\chi_e^2}
 {4\cosh^2((\omega_e+t\chi_e)/2)}\,dt
 =D(\operatorname{Ber}(p_e(1))\Vert
    \operatorname{Ber}(p_e(0))).
 }                                                        \tag{R34.3}
\]

Its closed form is

\[
 \boxed{
 K_e=\frac{\chi_e}{2}\tanh\frac{\omega_e+\chi_e}{2}
 -\log\frac{\cosh((\omega_e+\chi_e)/2)}
                {\cosh(\omega_e/2)}.
 }                                                        \tag{R34.4}
\]

This expression is unchanged when the temporary edge orientation is
reversed.  Thus (R34.3)--(R34.4) have no hidden sign convention.

Call an edge **crossing** when
`omega_e(omega_e+chi_e)<=0`.  On a noncrossing edge put

\[
 d_e=\min\{|\omega_e|,|\omega_e+\chi_e|\}.
\]

Since `(4 cosh^2(u/2))^(-1)<=e^(-|u|)`, (R34.3) gives the
useful rigorous split

\[
 K_e\le \frac{\chi_e^2}{2}e^{-d_e}\quad\hbox{off crossings},
 \qquad
 K_e\le\log(1+e^{|\omega_e|})\quad\hbox{on crossings}.   \tag{R34.5}
\]

The second inequality is the maximum possible KL from a Bernoulli law to the
fixed base `Ber(p_e(0))`.  The exact value (R34.4), rather than the upper
bound, should normally be retained.

## 2. The exact moving-edge-mass correction

For a fixed coordinate `j`, its edges partition the state space, so
`sum_e M_e(t)=1`.  Write expectation and covariance over this edge law as
`E_(M_(t,j))` and `Cov_(M_(t,j))`.  Since

\[
 \frac{d}{dt}\log M_e(t)
 =A'_e(t)-\mathbb E_{\mu_t}g,
\]

(R34.2) yields

\[
 \frac{d}{dt}\sum_{e\in E_j}M_e(t)k_e(t)
 =t\mathcal E_{t,j}(g)
 +\operatorname{Cov}_{M_{t,j}}(A'_e(t),k_e(t)).           \tag{R34.6}
\]

At `t=0` every `k_e` vanishes.  At `t=1`, the edge average is the
endpoint conditional relative entropy

\[
 C_j:=\sum_{e\in E_j}M_e(1)K_e
 =\mathbb E_{\mu_{1,-j}}
 D(\mu_1(D_j\mid D_{-j})\Vert
   \nu_\beta(D_j\mid D_{-j})).                            \tag{R34.7}
\]

Consequently

\[
 \boxed{
 \int_0^1t\mathcal E_{t,j}(g)\,dt
 =C_j-\int_0^1
 \operatorname{Cov}_{M_{t,j}}(A'_e(t),k_e(t))\,dt.
 }                                                        \tag{R34.8}
\]

Summing (R34.8) over the orientation coordinate and every projective vertex
coordinate is an exact formula for the full weighted energy in (10.961).
The covariance is the precise price of edge-context mass migrating during
the global interpolation.  Holding `M_e` fixed, or silently replacing it by
`M_e(1)`, deletes this term and is invalid.

A convenient sufficient upper bound following from (R34.8) is

\[
 \int_0^1t\mathcal E_t(g)\,dt
 \le \sum_j C_j+
 \int_0^1\sum_j
 [ -\operatorname{Cov}_{M_{t,j}}(A'_e(t),k_e(t)) ]_+\,dt. \tag{R34.9}
\]

This is sharper than independently taking suprema of `M_e(t)`: it preserves
the exact endpoint binary cost and isolates only the adverse part of the
mass migration.

## 3. Vertex resonance is posterior selector concentration

Use the canonical mixture from (10.938):

\[
 f(d)=\sum_Sq(S)h_S(d),\qquad
 \pi_d(S)=\frac{q(S)h_S(d)}{f(d)}.
\]

For a vertex `i`, define the posterior exclusion probability

\[
 r_i(d)=\pi_d\{S:i\notin S\}.
\]

If `y=T_i x`, every component with `i notin S` is unchanged because its
restriction does not see the flip.  Therefore

\[
 f(x)r_i(x)=\sum_{S:i\notin S}q(S)h_S(x)
            =\sum_{S:i\notin S}q(S)h_S(y)=f(y)r_i(y).
\]

Thus the full parent response has the exact signing-channel form

\[
 \boxed{
 \chi_i(x,y)=\log\frac{f(y)}{f(x)}
 =\log\frac{r_i(x)}{r_i(y)}.
 }                                                        \tag{R34.10}
\]

This is stronger than the one-sided moment-generating formula (10.939): it
uses the invariant selector subfamily specific to principal restriction.

The base odds remain

\[
 \omega_i(d)=-4\beta h_i(d).
\]

If `omega_i>0` but `omega_i+chi_i<=0`, then (R34.10) and `r_i(y)<=1`
give

\[
 r_i(x)\le e^{\chi_i}\le e^{-\omega_i}.
\]

If `omega_i<0` but `omega_i+chi_i>=0`, the symmetric conclusion is
`r_i(y)<=e^{-|omega_i|}`.  Hence

\[
 \boxed{
 \text{a vertex crossing forces posterior inclusion at the newly favored
 endpoint with probability at least }1-e^{-|\omega_i|}.
 }                                                        \tag{R34.11}
\]

There is also an exact mixture interpretation.  Under the endpoint joint law
`P(S,d)=mu_1(d)pi_d(S)`, condition on the edge context `D_-i` and put
`rho_e=P(i notin S|D_-i)`.  Given `i notin S`, the bit conditional is exactly
the base `Ber(p_e(0))`.  Thus, for some included-selector conditional `q_e`,

\[
 p_e(1)=\rho_ep_e(0)+(1-\rho_e)q_e,
\]

and convexity gives

\[
 K_e\le(1-\rho_e)D(\operatorname{Ber}(q_e)\Vert
                    \operatorname{Ber}(p_e(0)))
 \le(1-\rho_e)\log(1+e^{|\omega_e|}).                   \tag{R34.12}
\]

Moreover `sum_i E_(M_(1,i))(1-rho_e)=m` if all `n` vertex flips are counted.
This does not by itself have the target scale, but it identifies the exact
additional minimizer-specific tail needed: large base fields can contribute
at resonance only where the selector posterior has almost completely spent
its mass on selectors containing that vertex.

The orientation coordinate has no omitted-selector subfamily, so
(R34.10)--(R34.12) do not apply to it.  Its exact formulas remain

\[
 \omega_0(d)=-2\beta H_A(d),\qquad
 Z_{S,0}(d)=F_{\beta,S}(d[S])-F_{-\beta,S}(d[S]),
\]

with `chi_0=log E_(pi_d)e^(Z_(S,0))`.  It must be included in
`C_0`, the migration term, and the restoring drift below.

## 4. Exact restoring fit and its sharp scalar diagnostic

For the summed-coordinate heat-bath generator, the explicit drift is

\[
 \boxed{
 B_t(d)=-L_tg(d)
 =\sum_{j=0}^{n-1}
 \frac{\mu_t(T_jd)}{\mu_t(d)+\mu_t(T_jd)}
 [g(d)-g(T_jd)].
 }                                                        \tag{R34.13}
\]

This includes the orientation summand `j=0`.  Put

\[
 X_t=g-\mathbb E_{\mu_t}g,\quad
 V_t=\mathbb E X_t^2,\quad
 W_t=\mathbb E B_t^2,\quad
 E_t=\mathcal E_t(g)=\mathbb E X_tB_t.
\]

For arbitrary `kappa>0`, the restoring defect is exactly

\[
 R_t(\kappa)=W_t-2\kappa E_t+\kappa^2V_t.                \tag{R34.14}
\]

When `V_t,W_t,E_t>0`, define

\[
 \widehat\kappa_t=\sqrt{W_t/V_t},\qquad
 \rho_t=\frac{E_t}{\sqrt{V_tW_t}}\in(0,1].
\]

The choice `widehat kappa` minimizes the normalized defect
`R_t(kappa)/(kappa E_t)`, and direct substitution gives

\[
 \boxed{
 \frac{R_t(\widehat\kappa_t)}{\widehat\kappa_tE_t}
 =2(\rho_t^{-1}-1),
 \qquad
 C_{\rm scr}(t)=\frac{V_t}{E_t}
 =\frac1{\widehat\kappa_t\rho_t}.
 }                                                        \tag{R34.15}
\]

Thus a lower bound on the restoring rate `widehat kappa` and on the
alignment `rho` is the sharp scalar version of the endpoint-specific
comparison.  It is a diagnostic identity, not a proof that either lower
bound follows from exact minimality.

## 5. A fully quantified surviving sufficient lemma

Let `N_*=n^(1/2-2c)` be the parent target scale.  Uniformly over every
relevant order, exact minimizing signing, fixed-density target size, the
prescribed `beta=Theta(n^(-1/2+c))`, and `t in [0,1]`, it is enough to prove:

1. **Restoring comparison:** whenever `E_t(g)>0`, either (10.961) with constants
   `kappa_t>=kappa_0>0` and `D<infinity`, or the sharper scalar package

   \[
   \widehat\kappa_t\ge\kappa_0,
   \qquad \rho_t\ge\rho_0>0.                             \tag{R34.16}
   \]

2. **Endpoint resonant cost:** with the orientation edge included,

   \[
   \sum_{j,e}M_e(1)
   D(\operatorname{Ber}(p_e(1))\Vert
     \operatorname{Ber}(p_e(0)))=O(N_*).                 \tag{R34.17}
   \]

3. **No adverse mass migration above target scale:**

   \[
   \int_0^1\sum_j
   [-\operatorname{Cov}_{M_{t,j}}(A'_e(t),k_e(t))]_+dt
   =O(N_*).                                               \tag{R34.18}
   \]

Then (R34.8)--(R34.9) give
`int_0^1 t E_t(g)dt=O(N_*)`.  Under (R34.16), (R34.15) and (10.900) give

\[
 \operatorname{Ent}_{\nu_\beta}(f_\beta)
 \le(\kappa_0\rho_0)^{-1}
 \int_0^1t\mathcal E_t(g)dt=O(N_*).                     \tag{R34.19}
\]

Under the original defect formulation (10.961), the multiplier is exactly

\[
 \kappa_0^{-1}
 \left(\frac{\sqrt D+\sqrt{D+4}}2\right)^2.
\]

Thus the constants and the `t` weight match (10.900); no unweighted
substitution is hidden.  Equations (R34.5) and (R34.11) further reduce the
vertex part of (R34.17) to a nonresonant screened tail plus a weighted tail of
states with exponentially concentrated selector inclusion.  Establishing
that tail, (R34.18), and the orientation part requires actual
signing-minimizer structure.

## 6. Exhaustive finite numerical audit and obstructions

The checker reconstructs `f`, every posterior, all coordinate edges,
`B_t`, the endpoint conditional KL, and the covariance correction.  It
checks (R34.1), (R34.4), (R34.6)--(R34.10), the quadratic base odds, and
integration by parts to floating error below `2e-8` after quadrature.

The following decimals are **Numerical**, despite exhaustive finite-state
enumeration.

- On the exact `A_9` minimizer with `m=3`, `beta=0.5`, the maximum sampled
  `C_scr(t)` is about `4.894`, the minimum Rayleigh rate `E_t/V_t` is about
  `0.204`, and the maximum fixed-`kappa=1` defect ratio `R_t/E_t` is about
  `3.951`.  At `beta=2`, these become approximately `6.917`, `0.1446`, and
  `6.003`.  Hence the small `A_4/A_6` grid from Wave 33 does not support a
  small universal constant; an unspecified finite constant remains open.
- For `A_9,m=3,beta=2`, about `99.17%` of the integrated **vertex** energy is
  on crossing edges.
  A proof cannot treat resonance as a uniformly negligible exceptional set,
  even on known exact minimizers at finite order.  This fixed-order
  low-temperature wall does not match the project scaling by itself.
- For `A_9,m=7,beta=2`,

  \[
  \int_0^1t\mathcal E_tdt\approx0.05233614,
  \quad \sum_jC_j\approx0.05086819,
  \quad \int_0^1\sum_j\operatorname{Cov}\,dt
  \approx-0.00146794.
  \]

  Thus the migration covariance can have the adverse sign, and the shortcut
  `weighted energy <= endpoint conditional KL` is false on this finite exact
  minimizer.
- The earlier `A_6,beta=0.1,m=3` audit is reproduced: the orientation
  coordinate contributes about `34%` of the integrated energy.  In other
  examples orientation may be tiny, but it cannot be deleted algebraically.

## Frontier

The sharpest endpoint reduction found here is (R34.8): a sum of exact
resonant Bernoulli KL costs minus a single, explicitly identified
edge-context migration covariance.  The new project-specific content is
(R34.10)--(R34.11): a vertex resonance is exactly a posterior selector-
exclusion ratio and forces near-certain inclusion at the bias-reversing
endpoint.  This suggests a concrete minimizer-specific attack on weighted
tails of `-log r_i`, while (R34.17)--(R34.18) state precisely what is needed
to close the parent energy.  The orientation coordinate and the restoring
angle/rate remain independent inputs.  `A_9` rules out rare-resonance,
favorable-migration-sign, and small-fixed-constant shortcuts, but supplies no
asymptotic falsifier of the surviving lemma.
