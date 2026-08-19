# Zero-regret switching endpoints have an orthogonal bridge covariance

Status: **proved exact endpoint identity and scalable method-class
obstruction for actual optimizing children**.  Biased vertex switching does
remove the random-competitor regrets from the Wishart--Wigner endpoint.
However, its complete covariance lies in the sector which is even under the
two child global-spin flips, whereas the `mn` bridge features lie in the
sector which is odd under both flips.  Consequently no sign-definite
covariance comparison can turn this zero-regret endpoint into a sublinear
composition bound without adding a correction of rank at least `mn` and
total variance at least `t^2mn=Theta_beta(N)` on comparable splits.

This is a no-go for the raw zero-regret Guerra/Slepian architecture, not for
a genuinely signed interpolation.  In particular, the result does not rule
out a law-specific self-saturation estimate along the actual interpolating
Gibbs measures.  Section 5 does prove that merely retaining the replica sign
and then taking a supremum over replica laws still has a linear floor, even
when every favorable local-field covariance is kept.

## 1. Setup and the direct defect implication

Let

```math
\phi_A(u)=\log\mathbb E_x\cosh(uH_A(x)),
\qquad N=m+n,
\qquad t={\beta\over\sqrt N},
```

and put `s_m=beta/sqrt(m)`, `s_n=beta/sqrt(n)`.  Let `A,D` be exact
minimizers defining `P_m(beta),P_n(beta)`, and let

```math
\mu_A(\tau,x)
={2^{-m-1}e^{t\tau H_A(x)}\over e^{\phi_A(t)}},
\qquad
\mu_D(\sigma,y)
={2^{-n-1}e^{t\sigma H_D(y)}\over e^{\phi_D(t)}}.
\tag{1.1}
```

For an orientation `epsilon`, the exact bridge endpoint is

```math
\mathcal F_B^\epsilon
=\mathbb E_{B}\log\left[
2\mathbb E_{\mu_A\mu_D}{\bf1}_{\{\tau\sigma=\epsilon\}}
e^{t\tau x^{\mathsf T}By}\right],
\qquad
\mathcal F_B={\mathcal F_B^++\mathcal F_B^-\over2}.
\tag{1.2}
```

The exact selected-child construction identity is

```math
E_{m,n}(\beta)\le
\mathcal F_B-\Delta_A-\Delta_D,
\qquad
\Delta_A=\phi_A(s_m)-\phi_A(t),
\quad
\Delta_D=\phi_D(s_n)-\phi_D(t).
\tag{1.3}
```

Set

```math
\lambda_m={t\over s_m}=\sqrt{m/N},
\qquad
\lambda_n={t\over s_n}=\sqrt{n/N}.
\tag{1.4}
```

Let `z_1,...,z_m` be independent signs of mean `sqrt(lambda_m)` and
let `w_1,...,w_n` be independent signs of mean `sqrt(lambda_n)`.
Thus `E z_i z_j=lambda_m` and `E w_a w_b=lambda_n` off the diagonal.
Write

```math
A^z=\operatorname{diag}(z)A\operatorname{diag}(z),
\qquad
D^w=\operatorname{diag}(w)D\operatorname{diag}(w),
\tag{1.5}
```

and define the switching endpoint

```math
\begin{aligned}
\mathcal F_{\rm sw}^\epsilon
=\mathbb E_{z,w}\log\bigg[
2\mathbb E_{\mu_A\mu_D}{\bf1}_{\{\tau\sigma=\epsilon\}}
\exp\{&\tau(s_mH_{A^z}(x)-tH_A(x))\\
       &+\sigma(s_nH_{D^w}(y)-tH_D(y))\}\bigg],
\end{aligned}
\tag{1.6}
```

with `F_sw` the average of the two orientations.  Every realization in
(1.5) is gauge-equivalent to its child.  A change of variables
`x -> z*x`, `y -> w*y` in (1.6) therefore gives the exact, realization-wise
formula

```math
\mathcal F_{\rm sw}^\epsilon
=\Delta_A+\Delta_D
+\log\{1+\epsilon u_A(s_m)u_D(s_n)\},
\tag{1.7}
```

where

```math
u_A(s)={\mathbb E_x\sinh(sH_A(x))\over
              \mathbb E_x\cosh(sH_A(x))}
```

and similarly for `D`.  Hence

```math
\boxed{
\mathcal F_{\rm sw}
=\Delta_A+\Delta_D
+{1\over2}\log\{1-u_A(s_m)^2u_D(s_n)^2\}
\le\Delta_A+\Delta_D.}
\tag{1.8}
```

This endpoint has literally zero child regret.  In particular, (1.3) and
(1.8) give the immediate quantitative arrow

```math
\boxed{
\mathcal F_B\le\mathcal F_{\rm sw}+\omega_N
\quad\Longrightarrow\quad
E_{m,n}(\beta)\le\omega_N.}
\tag{1.9}
```

Thus comparing these two endpoints to power-saving accuracy would prove the
requested composition bound directly; no new scalar surrogate is being
introduced.

## 2. Exact covariance of the switching perturbation

For the left child define

```math
V_A^z(x)
=s_m\sum_{i<j}a_{ij}(z_iz_j-\lambda_m)x_ix_j,
\qquad
h_A(x)_i=x_i(Ax)_i.
\tag{2.1}
```

For two replicas put `X=<x^1,x^2>`.  Then

```math
\boxed{
\mathbb E_zV_A^z(x^1)V_A^z(x^2)
=s_m^2\left[
(1-\lambda_m)^2{X^2-m\over2}
+\lambda_m(1-\lambda_m)
 \langle h_A(x^1),h_A(x^2)\rangle
\right].}
\tag{2.2}
```

Indeed, for `e={i,j}` and `f={k,l}`,

```math
\operatorname{Cov}(z_iz_j,z_kz_l)
=\begin{cases}
1-\lambda_m^2,&e=f,\\
\lambda_m(1-\lambda_m),&e\ne f\text{ and }|e\cap f|=1,\\
0,&e\cap f=\varnothing.
\end{cases}
\tag{2.3}
```

Moreover,

```math
\sum_i h_A(x^1)_ih_A(x^2)_i
=2\sum_e x_e^1x_e^2
+\sum_{\substack{e\ne f\\|e\cap f|=1}}
a_ea_f x_e^1x_f^2,
\tag{2.4}
```

where `x_{ij}^r=x_i^rx_j^r`, while
`sum_e x_e^1x_e^2=(X^2-m)/2`.  Substitution of (2.3)--(2.4) proves
(2.2).  The right child has the identical formula with `n,D,Y`.

After multiplication by the common replica sign in one orientation sector,
the full switching covariance is

```math
\begin{aligned}
C_{\rm sw}(1,2)=\tau^1\tau^2\{&s_m^2(1-\lambda_m)^2{X^2-m\over2}
+s_m^2\lambda_m(1-\lambda_m)
 \langle h_A(x^1),h_A(x^2)\rangle\\
&+s_n^2(1-\lambda_n)^2{Y^2-n\over2}
+s_n^2\lambda_n(1-\lambda_n)
 \langle h_D(y^1),h_D(y^2)\rangle\}.
\end{aligned}
\tag{2.5}
```

By contrast, the covariance-matched Gaussian bridge has

```math
\boxed{C_B(1,2)=t^2\tau^1\tau^2XY.}
\tag{2.6}
```

Formula (2.5) also displays a concrete Gaussian realization: independent
child-edge Gaussians with coefficient `s_k(1-lambda_k)` and independent
Gaussian coordinates multiplying the local-field vectors with coefficient
`s_k sqrt(lambda_k(1-lambda_k))`.

The exact zero-regret statement (1.8) concerns the Bernoulli switching
chaos, whereas (2.5) is its covariance-matched Gaussian process.  Passing
between those two free energies would require a separate universality
estimate.  The obstruction below occurs even after granting that step for
free.

## 3. The parity-sector obstruction

Let the state space in one orientation sector be

```math
\Omega_\epsilon
=\{(\tau,x,y):\tau\in\{\pm1\},\ x\in\{\pm1\}^m,
 y\in\{\pm1\}^n\},
\tag{3.1}
```

with `sigma=epsilon*tau`.  The Gaussian feature functions generating
(2.5) are

```math
\tau x_ix_j,
\quad \tau h_A(x)_i,
\quad \tau y_ay_b,
\quad \tau h_D(y)_a.
\tag{3.2}
```

Every function in (3.2) is invariant under `x -> -x` and under
`y -> -y` separately.  The bridge features

```math
b_{ia}(\tau,x,y)=\tau x_i y_a                              \tag{3.3}
```

change sign under either global flip.  Thus the switching and bridge
covariances occupy orthogonal character sectors.  This is not a matter of
the signs or spectra of `A,D`.

The obstruction is already visible in the canonical metrics.  For any
state `q=(tau,x,y)` and `q'=(tau,-x,y)`, the switching process takes the
same value at `q,q'`, whereas the bridge process changes sign.  Therefore

```math
\boxed{
d_{\rm sw}(q,q')^2=0,
\qquad
d_B(q,q')^2
=\mathbb E\{G_B(q)-G_B(q')\}^2
=4t^2mn={4\beta^2mn\over N}.}
\tag{3.4}
```

On every split with `m,n>=kappa N`, the missing squared distance in (3.4)
is at least `4 beta^2 kappa^2 N`.

More generally, for every centered endpoint process `G_even` which is
separately invariant under the two child global flips,

```math
\boxed{
\sup_{q,q'\in\Omega_\epsilon}
\{d_B(q,q')^2-d_{\rm even}(q,q')^2\}_+
\ge {4\beta^2mn\over N}.}
\tag{3.4a}
```

Thus even granting arbitrary optimizer-dependent covariance inside the
separately-even sector, a uniform canonical-metric comparison has a linear
remainder on comparable splits.

There is also an all-directions statement.  Equip functions on
`Omega_epsilon` with the uniform inner product.  The `mn` functions in
(3.3) are orthonormal, and their span is orthogonal to every separately-even
feature, including all of (3.2).  Hence, if a positive-semidefinite
correction kernel `K` is required to make

```math
C_{\rm sw}+K-C_B\succeq0,                              \tag{3.5}
```

then compression of (3.5) to `span{b_ia}` gives

```math
K\big|_{\operatorname{span}\{b_{ia}\}}
\succeq t^2I_{mn}.
```

Consequently

```math
\boxed{
\operatorname{rank}K\ge mn,
\qquad
\mathbb E_{q\sim U}K(q,q)\ge t^2mn
={\beta^2mn\over N}.}
\tag{3.6}
```

The same conclusion holds if one augments (2.5) by any number of processes
which remain separately even in the two child spins.

## 4. Quantitative no-go and exact scope

Equations (1.8)--(1.9) show why the switching endpoint was attractive: a
comparison error `O(N^(1-delta))` would immediately give the desired defect.
Equations (3.4)--(3.6) show that the standard sign-definite route to such a
comparison cannot work:

```text
zero-regret child switching endpoint
    + separately-even covariance corrections
    + positive-semidefinite Guerra/Slepian domination
```

requires `mn` new odd bridge directions and a variance budget
`beta^2mn/N=Theta_beta(N)` on every comparable split.  In particular, a
comparison which pays a positive-semidefinite correction by its universal
annealed variance bound retains a linear certificate and cannot establish
(1.9) with a power saving.

The obstruction applies to arbitrary children and therefore, in particular,
to every selected pair of exact own-scale pressure minimizers, at every
order and every split.  It is stronger than saying that the overlap
quadratic form in (2.5)--(2.6) is indefinite: it identifies an entire
`mn`-dimensional parity sector on which the zero-regret endpoint has no
covariance at all.

Two escapes remain outside the theorem.

1. A signed interpolation may keep the replica factor `tau^1 tau^2` and
   obtain cancellation without positive-semidefinite covariance ordering.
2. One may add genuinely joint, odd-in-each-child processes and control
   their quenched free energy more sharply than by total variance.  Such a
   correction is already bridge-like: (3.6) says that it must supply all
   `mn` odd directions, not a scalar gauge mode or a bounded collection of
   child statistics.

Thus biased switching removes the optimizer regrets exactly, but gauge
equivalence also removes precisely the covariance sector which the bridge
occupies.  The raw covariance-matched zero-regret interpolation is therefore
not a sublinear recurrence architecture.

## 5. Signed interpolation: exact derivative and a measure-uniform no-go

The replica sign does permit a non-positive-semidefinite comparison, but it
does not give a uniform closure over Gibbs laws.  This can be stated without
discarding any favorable child term.

Let `G_B` be the Gaussian process with covariance (2.6), let `G_sw` be the
Gaussian process with covariance (2.5), and interpolate

```math
G_u=\sqrt u\,G_B+\sqrt{1-u}\,G_{\rm sw}.
\tag{5.1}
```

Use the finite base measure
`2 1_{tau sigma=epsilon} mu_A mu_D` in either orientation sector.  Write
`Var_u` for Gibbs variance followed by disorder expectation, and set

```math
e_{ij}=\tau x_ix_j,
\quad f_i=\tau h_A(x)_i,
\quad
\widetilde e_{ab}=\tau y_ay_b,
\quad \widetilde f_a=\tau h_D(y)_a,
\tag{5.2}
```

where the harmless fixed orientation sign has been absorbed on the right.
Gaussian integration by parts gives the exact identity

```math
\boxed{
\begin{aligned}
F_B^G-F_{\rm sw}^G
={1\over2}\int_0^1\bigg[&t^2\sum_{i,a}\operatorname{Var}_u(\tau x_iy_a)\\
&-a_m\sum_{i<j}\operatorname{Var}_u(e_{ij})
-b_m\sum_i\operatorname{Var}_u(f_i)\\
&-a_n\sum_{a<b}\operatorname{Var}_u(\widetilde e_{ab})
-b_n\sum_a\operatorname{Var}_u(\widetilde f_a)
\bigg]\,du,
\end{aligned}}
\tag{5.3}
```

with

```math
a_m=s_m^2(1-\lambda_m)^2,
\qquad b_m=s_m^2\lambda_m(1-\lambda_m),
\tag{5.4}
```

and the analogous right-child coefficients.  Thus all local-field and
child-edge terms have the favorable sign.  The only positive term is the
variance in the odd bridge sector exposed in Section 3.

Formula (5.3) has a direct defect implication.  If

```math
|\mathcal F_B-F_B^G|
+|F_{\rm sw}^G-\mathcal F_{\rm sw}|\le\eta_N
\tag{5.5}
```

after averaging the two orientations, and the right side of (5.3) is at
most `omega_N`, then (1.9) gives

```math
\boxed{E_{m,n}(\beta)\le\eta_N+\omega_N.}             \tag{5.6}
```

The bridge part of `eta_N` is `O_beta(sqrt(N))` by the usual
coordinatewise Lindeberg replacement.  No sublinear bound for the
switching-chaos part is asserted here.

There is a sharp obstruction to closing (5.3) by an inequality valid for
all replica laws.  For an arbitrary state `q=(tau,x,y)`, let `nu_q` be the
uniform law on

```math
q=(\tau,x,y),
\qquad q'=(\tau,-x,y).                                \tag{5.7}
```

Every child-edge and local-field feature in (5.2) has zero variance under
`nu_q`, while every bridge feature `tau x_i y_a` has variance one.  Hence
the variance functional in the integrand of (5.3), evaluated at `nu_q`, is
exactly

```math
\boxed{
\mathscr D(\nu_q)={t^2mn\over2}
={\beta^2mn\over2N}.}                                \tag{5.8}
```

Consequently, for every split and every pair of children,

```math
\boxed{
\sup_{\nu\in\mathcal P(\Omega_\epsilon)}\mathscr D(\nu)
\ge {\beta^2mn\over2N}.}                             \tag{5.9}
```

On comparable splits this is `Theta_beta(N)`.  Therefore any signed
interpolation argument which forgets the actual interpolating Gibbs law and
upper-bounds (5.3) by a measure-uniform covariance or overlap inequality
still has a linear floor.  The replica sign does not help this witness,
because the two states in (5.7) use the same `tau`.

The obstruction is not an artifact of suppressing all child variation.
Fix `tau` and `y`, and let `x` be uniform on its entire Boolean cube.  Call
this law `nu_left`.  Then

```math
\sum_{i,a}\operatorname{Var}_{\nu_{\rm left}}(\tau x_iy_a)=mn,
\qquad
\sum_{i<j}\operatorname{Var}_{\nu_{\rm left}}(e_{ij})=K_m,
\tag{5.10}
```

and, for every signing `A`,

```math
\sum_i\operatorname{Var}_{\nu_{\rm left}}(f_i)
=\mathbb E_x\|Ax\|_2^2=m(m-1)=2K_m.                 \tag{5.11}
```

All right-child variances vanish.  The switching coefficients satisfy the
exact cancellation

```math
a_m+2b_m=s_m^2(1-\lambda_m^2)=s_m^2-t^2.             \tag{5.12}
```

Therefore the full signed derivative functional, with every favorable
local-field term retained, is

```math
\boxed{
\mathscr D(\nu_{\rm left})
={1\over2}\{t^2mn-(s_m^2-t^2)K_m\}
={\beta^2n(m+1)\over4N}.}                            \tag{5.13}
```

This is again `Theta_beta(N)` on every comparable split (and is
`beta^2(N+2)/16` at an equal split).  Thus neither the replica sign nor the
complete local-field covariance supplies a measure-uniform cancellation.
The obstruction is independent of `A,D` and hence includes the actual
optimizing children.

The only surviving signed route is law-specific self-saturation: one must
prove that the actual Gibbs measures generated by (5.1) suppress the odd
bridge variance strongly enough, in integral over `u`, to offset (5.8),
while retaining the favorable child variances.  Equation (5.6) gives the
exact quantitative target.  A supremum-over-overlaps closure, even one
retaining all the local-field terms in (5.3), cannot provide it.

## 6. Optimizer edge-flip stationarity alone still has a linear model

One might hope to control the actual Gibbs laws in (5.3) using the child
optimizer inequalities.  At inverse temperature `s`, flipping one child
edge gives the necessary stationarity condition

```math
a_{ij}\mathbb E_\nu(\tau x_ix_j)\le\tanh s.           \tag{6.1}
```

Indeed, the exact flip ratio is

```math
\mathbb E_\nu e^{-2s\tau a_{ij}x_ix_j}
=\cosh(2s)-a_{ij}\mathbb E_\nu(\tau x_ix_j)\sinh(2s),
\tag{6.2}
```

and a minimizing signing requires (6.2) to be at least one.  The exact
own-temperature Gibbs law also has

```math
\sum_{i<j}a_{ij}\mathbb E_\nu(\tau x_ix_j)
=\phi_A'(s)\ge0.                                    \tag{6.3}
```

Even granting (6.1) and (6.3) for both interpolating child marginals does
not imply a sublinear bound in (5.3).  Here is an explicit
arbitrarily-large-order replica-law countermodel.

Take an even equal child order `r`, put `A=D=J-I`, fix the positive
orientation, and let `tau=sigma=-1`.  Independently choose `x,y` uniformly
from the balanced slice

```math
\{x\in\{\pm1\}^r:\ \sum_i x_i=0\}.                  \tag{6.4}
```

For distinct coordinates under (6.4),

```math
\mathbb E x_ix_j=-{1\over r-1}.
\tag{6.5}
```

Thus every aligned edge mean in (6.1) equals `1/(r-1)`, and (6.3) equals
`r/2>0`.  Consequently all one-edge stationarity inequalities for both
children hold whenever

```math
{1\over r-1}\le\tanh(\beta/\sqrt r),                 \tag{6.6}
```

which is true for every fixed `beta>0` and all sufficiently large even
`r` (the ratio of the right side to the left side is asymptotic to
`beta sqrt(r)`).

The local fields vanish as random quantities in the derivative.  In fact,
for every balanced `x`,

```math
h_A(x)_i=x_i\sum_{j\ne i}x_j=-1,
\qquad \tau h_A(x)_i=1,                              \tag{6.7}
```

and similarly on the right.  Hence all local-field variances in (5.3) are
zero.  The edge variances are

```math
\operatorname{Var}(\tau x_ix_j)
=1-{1\over(r-1)^2},                                  \tag{6.8}
```

whereas every bridge feature has mean zero and variance one.  At the equal
split, `lambda=2^{-1/2}` and

```math
a=(3-2\sqrt2)t^2.
\tag{6.9}
```

Substitution in the complete signed functional from (5.3) gives the exact
finite-order value

```math
\boxed{
\begin{aligned}
\mathscr D_r^{\rm bal}
&={1\over2}\left[
t^2r^2-2aK_r\left(1-{1\over(r-1)^2}\right)\right]\\
&={t^2r^2\over2}\left[
1-(3-2\sqrt2){r-2\over r-1}\right].
\end{aligned}}                                      \tag{6.10}
```

Since `N=2r` and `t^2=beta^2/N`,

```math
\boxed{
{\mathscr D_r^{\rm bal}\over N}
\longrightarrow {\sqrt2-1\over4}\,\beta^2>0.}      \tag{6.11}
```

This is a countermodel to a **stationarity-only closure**, not a claim that
`J-I` is a pressure optimizer or that (6.4) is the actual interpolating
Gibbs law.  Its quantifier is the relevant one: every inequality derived
solely from the necessary optimizer constraints (6.1), the radial sign
(6.3), Boolean support, and the exact local-field definitions must also
hold for (6.4), but (6.10) forces its upper bound to be linear.  Therefore
optimizer edge-flip stationarity cannot by itself turn (5.6) into a
power-saving defect.  Any successful law-specific argument must use a
stronger property of the actual interpolating Gibbs law than all
coordinatewise flip inequalities and nonnegative radial energy.
