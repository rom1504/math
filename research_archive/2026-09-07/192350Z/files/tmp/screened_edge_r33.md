# Wave 33 route 3: edge-symmetrized screened response

## Status

- **Verified identity:** on every parent coordinate edge, the response has an
  unoriented formula involving both posterior KL directions.  The exact
  heat-bath cost is the square of this full signed combination, multiplied by
  a conductance which depends on the base Gibbs odds and the response.
- **Verified identity:** the global orientation edge has a different residual
  from a vertex edge.  Its base log-odds is `-2 beta H_A(d)` and its selector
  residual is the difference between the `+beta` and `-beta` outside free
  energies.  It cannot in general be discarded.
- **Verified sufficient lemma:** a one-function heat-bath drift-defect bound
  controls the endpoint comparison coefficient without asserting a generic
  Poincare inequality.  Together with a target-scale bound on the exact
  resonant screened energy, it proves the parent entropy target.
- **Finite mechanism walls:** `C_scr <= 1` and pointwise restoring drift both
  fail on exact small minimizers.  Conversely, an arbitrary positive endpoint
  likelihood can have `C_scr -> infinity`, even with duplicate selectors and
  zero posterior KL.  Thus a bounded one-function comparison still needs
  harmonic/minimizer structure.
- **Open at project scale:** neither the drift-defect bound nor the resonant
  energy bound is proved uniformly for exact minimizers.  The separate
  adjacent-selector Johnson--Hellinger estimate is untouched.

The exact finite checks are in `tmp/screened_edge_r33.py`.

## 1. Both posterior KL directions give the exact unoriented response

Fix an unordered coordinate edge `e={x,y}`, temporarily orient it from `x`
to `y`, and write

\[
 Z_e(S)=\log\frac{h_S(y)}{h_S(x)},\qquad
 \chi_e=\log\frac{f_\beta(y)}{f_\beta(x)}.
\]

Let `pi_x,pi_y` be the selector posteriors and put

\[
 K_x=D(\pi_x\Vert\pi_y),\quad K_y=D(\pi_y\Vert\pi_x),\quad
 a_x=\mathbb E_{\pi_x}Z_e,\quad a_y=\mathbb E_{\pi_y}Z_e.
\]

The exponential tilt `pi_y=pi_x exp(Z_e-chi_e)` gives, in both directions,

\[
 \boxed{
 a_x=\chi_e-K_x,\qquad
 a_y=\chi_e+K_y.
 }
\tag{R33.1}
\]

Consequently

\[
 \boxed{
 \chi_e=\frac12\{a_x+a_y+K_x-K_y\},\qquad
 a_y-a_x=K_x+K_y.
 }
\tag{R33.2}
\]

The second identity is the posterior Jeffreys divergence.  Reversing the edge
interchanges `x,y`, replaces `Z_e,chi_e` by their negatives, and changes the
first bracket in (R33.2) by a sign.  Thus its square is canonically attached
to the unordered edge.

This is sharper than separately using `chi=a_x+K_x`: the difference
`K_x-K_y` in (R33.2) can cancel the symmetrized drift `a_x+a_y`.  Triangle or
separate-square estimates erase precisely the cancellation seen in the
`A_6` audit from Wave 32.

## 2. Base odds, conductance, and the orientation coordinate

Let

\[
 \omega_e=\log\frac{\nu_\beta(y)}{\nu_\beta(x)},\qquad
 u_{t,e}=\log\frac{\mu_t(y)}{\mu_t(x)}
          =\omega_e+t\chi_e,qquad
 M_{t,e}=\mu_t(x)+\mu_t(y).
\]

The conditional probability of jumping from `x` to `y` in the coordinate
heat bath is `p_{t,e}(x)=1/(1+exp(-u_(t,e)))`.  The unordered conductance is

\[
 \boxed{
 c_t(e)=\frac{\mu_t(x)\mu_t(y)}{\mu_t(x)+\mu_t(y)}
 =\frac{M_{t,e}}{4\cosh^2(u_{t,e}/2)}.
 }
\tag{R33.3}
\]

Combining (R33.2)--(R33.3) gives the exact edge cost

\[
 \boxed{
 c_t(e)\chi_e^2
 =\frac{M_{t,e}}{16\cosh^2((\omega_e+t\chi_e)/2)}
   \{a_x+a_y+K_x-K_y\}^2.
 }
\tag{R33.4}
\]

Thus only a large response near the resonance
`omega_e+t chi_e approximately 0` is expensive.  Large aligned base and
endpoint slopes are automatically suppressed by the actual `mu_t`
conductance.  A useful estimate should target (R33.4) as a whole, rather than
bound the drift, two KL terms, response, and conductance independently.

For a projective vertex flip `T_i`, with

\[
 h_i(d)=\sum_{k\ne i}A_{ik}d_{ik},
\]

the base odds and selector residual are

\[
 \boxed{
 \omega_i(d)=-4\beta h_i(d),\qquad
 Z_{S,i}(d)=2\mathbf1_{\{i\in S\}}(b-\beta a).
 }
\tag{R33.5}
\]

The coordinate system also has the global orientation involution `T_0d=-d`.
Since `H_A(-d)=-H_A(d)` and negating every oriented extension negates its
outside energy,

\[
 F_{\beta,S}(-d[S])=F_{-\beta,S}(d[S]).
\]

Therefore the distinct orientation formulas are

\[
 \boxed{
 \omega_0(d)=-2\beta H_A(d),\qquad
 Z_{S,0}(d)=F_{\beta,S}(d[S])-F_{-\beta,S}(d[S]).
 }
\tag{R33.6}
\]

Equations (R33.1)--(R33.4) apply unchanged to this edge, but the cavity
formula in (R33.5) does not.

## 3. An endpoint-only comparison certificate

Put `g=log f_beta` and define the coordinate heat-bath generator at time `t`
by

\[
 (L_tg)(d)=\sum_{j=0}^{n-1}
 \frac{\mu_t(T_jd)}{\mu_t(d)+\mu_t(T_jd)}\,
 \chi_j(d)
 =\sum_j\frac{\chi_j(d)}{1+e^{-[\omega_j(d)+t\chi_j(d)]}}.
\tag{R33.7}
\]

Let `gbar_t=E_(mu_t)g`, `X_t=g-gbar_t`, and `B_t=-L_tg`.  Reversibility gives
the exact one-function integration-by-parts identity

\[
 \boxed{
 \mathcal E_t(g)=\sum_{j,e\in E_j}c_t(e)\chi_j(e)^2
 =\mathbb E_{\mu_t}X_tB_t.
 }
\tag{R33.8}
\]

When `E_t(g)>0`, this suggests a concrete comparison weaker than a Poincare
theorem.  Choose any scalar `kappa_t>0` and define the restoring-drift defect

\[
 R_t=\mathbb E_{\mu_t}(B_t-\kappa_tX_t)^2,qquad
 \delta_t=\frac{R_t}{\kappa_t\mathcal E_t(g)}.
\]

Writing `B_t=kappa_t X_t+r_t` in (R33.8), followed by one Cauchy--Schwarz
inequality, proves

\[
 \boxed{
 C_{\rm scr}(t)=\frac{\operatorname{Var}_{\mu_t}(g)}{\mathcal E_t(g)}
 \le\frac1{\kappa_t}
 \left(\frac{\sqrt{\delta_t}+\sqrt{\delta_t+4}}2\right)^2.
 }
\tag{R33.9}
\]

Indeed, if `Y^2=kappa_t Var(g)/E_t(g)`, then
`Y^2 <= 1+sqrt(delta_t)Y`.  This proves (R33.9), including its constants.
It is a certificate only for the actual endpoint `g`; it makes no assertion
about any other function or the spectral gap of `mu_t`.  If `E_t(g)=0`, the
connected coordinate cube forces `g` to be constant, so its variance and
parent contribution are zero and this ratio argument is unnecessary.

## 4. Exact project-scale sufficient target

At the fixed-density project temperature, it is enough to find constants
`kappa_0>0,D<infinity`, uniformly in `n`, every relevant exact minimizing
signing and target pair, every fixed-density deletion size in the project
window, and `t in [0,1]` (at the prescribed project-scale `beta`), such that
some `kappa_t>=kappa_0` obeys

\[
 \boxed{
 \mathbb E_{\mu_t}(B_t-\kappa_tX_t)^2
 \le D\kappa_t\mathcal E_t(g),
 }
\tag{R33.10}
\]

and simultaneously prove the resonant screened-energy estimate, with its
implicit constant uniform over the same data,

\[
 \boxed{
 \int_0^1t\sum_{j,e\in E_j}
 \frac{M_{t,e}}{16\cosh^2((\omega_e+t\chi_e)/2)}
 \{a_x+a_y+K_x-K_y\}^2dt
 =O(n^{1/2-2c}).
 }
\tag{R33.11}
\]

Equations (10.900), (R33.4), and (R33.9) then give the desired parent entropy
bound with the fixed multiplier

\[
 \kappa_0^{-1}
 \left(\frac{\sqrt D+\sqrt{D+4}}2\right)^2.
\]

This replaces the generic `C_HB(mu_t)` in (10.940) by an endpoint-specific
drift certificate while preserving the signed drift--KL cancellation and the
base-Gibbs conductance screening.  Both (R33.10) and (R33.11) are **open
minimizer-specific lemmas**, not consequences of the edge algebra.  They
address only the parent term; selector Hellinger is still separately needed.

## 5. Finite audit and scoped walls

Exact enumeration of `A_4` and `A_6`, for
`beta in {0.1,0.5,1,2,4}`, every available `m`, and
`t in {0,1/2,1}`, verifies (R33.1)--(R33.8) with maximum error below
`2.2e-14`.  The orientation/free-energy identity (R33.6) has error below
`4.5e-16`.

Three proposed shortcuts fail:

1. **`C_scr<=1` is false.**  On `A_4`, `beta=1/2,m=3,t=0`, exact enumeration
   gives `C_scr=1.398592860...`.
2. **Pointwise restoring drift is too strong.**  The condition
   `X_t(d)B_t(d)>=kappa X_t(d)^2` fails at three `A_4` states in the preceding
   case.  It also fails at twelve states for `A_6,beta=1/2,m=4,t=0`, where
   `min B_t/X_t=-9.29915983...`, even though `C_scr=0.371941426...`.
3. **The orientation edge cannot be deleted algebraically.**  For
   `A_6,beta=0.1,m=3,t=0`, it contributes about `34%` of the exact energy.

The L2 defect (R33.10) is substantially less brittle on these finite cases.
With `kappa_t=1`, the tested grid has

\[
 0.31719<C_{\rm scr}<1.39860,qquad
 R_t/\mathcal E_t(g)\le1.6123.
\]

This is numerical finite evidence only; it does not establish constants at
project scale.

Finally, no bounded one-function comparison follows from likelihood algebra.
Take the uniform base law on a two-dimensional cube and, up to normalization,

\[
 f=(1,\varepsilon,\varepsilon,e)
\]

at its four vertices.  At `t=1`, the variance of `log f` stays bounded away
from zero while its edge energy is
`Theta(epsilon log^2(1/epsilon))`; hence
`C_scr=Theta(1/(epsilon log^2(1/epsilon))) -> infinity`.  The checked values
for `epsilon=10^(-2),10^(-4),10^(-6)` are respectively
`1.2404,19.8558,891.2523`.  Duplicating this same likelihood across selectors
makes both posterior KL directions identically zero, so posterior algebra
does not remove the bottleneck.  This is an abstract endpoint wall, not a
harmonic endpoint of a complete-signing minimizer.

## Frontier

Edge symmetrization makes the surviving harmonic target more explicit.  The
dangerous quantity is not raw cavity response or posterior susceptibility;
it is the resonant, conductance-weighted square in (R33.4), with both KL
directions kept inside the signed bracket.  A genuinely endpoint-specific
replacement for generic parent Poincare is the drift-defect lemma (R33.10).
Finite minimizers support its L2 form but disprove the tempting unit-constant
and pointwise versions.  The route remains open until exact global minimality
or another signing-specific mechanism proves both (R33.10) and (R33.11) at
the target scale.
