# Wave 31 route 3: exact harmonic composition and the centered-accumulation wall

## Status

- **Verified:** the all-deletion harmonic cost is a reverse conditional KL.
  It has an exact associative block composition and an exact one-spin
  renormalized-cavity recursion.
- **Verified:** the endpoint soft minimum has an exact selector-posterior
  response formula.  A vertex-flip response is a soft average of the
  renormalized-field screening residual, not of the raw row field.
- **Verified wall:** entropy-chain, first-moment, and uncentered block
  accumulation bounds pay nonnegative reverse-KL common modes which the
  normalized endpoint likelihood cancels.  This loss is already drastic for
  the exact minimizer `A_6` at fixed density.
- **Not proved:** no fixed-density `O(n^(1/2-2c))` parent bound follows from
  exact minimality.  Doob gives a sharp centered criterion, while ordinary
  Efron--Stein needs a parent-specific dependence theorem and the available
  flip bound is off scale.

Every identity and finite value below is checked by
`tmp/harmonic_composition_r31.py`.

## 1. Exact one-spin integration

Write `V=[n]`, let

\[
E_A(d)=\langle A,d\rangle,
\qquad c_S(d)=\langle A[S],d[S]\rangle,
\]

and, for `y=d[S]` and `k=n-|S|`, normalize the outside partition sum as

\[
G^V_S(y)=\frac{K_{\beta,S}(y)}{2^k},
\qquad F^V_S(y)=\log G^V_S(y).
\tag{R31.1}
\]

Thus `F^V_V=0`.  Take `R=S union {i}`.  Let `d^+` and `d^-` be the two
oriented-cut extensions of `y` to `R`, chosen so that `d^-` is obtained by
flipping vertex `i`, and put

\[
a_{S,i}(d^+)=c_R(d^+)-c_S(y)
=2\sum_{j\in S}A_{ij}d^+_{ij}.
\tag{R31.2}
\]

The other increment is `-a`.  Directly partitioning the outside sum gives

\[
\boxed{
G^V_S(y)=\frac12\left(
 e^{\beta a}G^V_R(d^+)+e^{-\beta a}G^V_R(d^-)
\right).
}
\tag{R31.3}
\]

Put `F^+=F^V_R(d^+)`, `F^-=F^V_R(d^-)` and

\[
b^V_{S,i}(y)=\beta a+\frac{F^+-F^-}{2}.
\tag{R31.4}
\]

Then the exact nonlinear recursion is

\[
\boxed{
F^V_S(y)=\frac{F^++F^-}{2}+\log\cosh b^V_{S,i}(y).
}
\tag{R31.5}
\]

The parent conditional marginal of the new spin is

\[
\nu_\beta(d^\pm\mid y)
=\frac{e^{\pm b}}{2\cosh b}.
\tag{R31.6}
\]

Thus `b` is the exact renormalized cavity half-log-odds after all remaining
outside spins have been integrated.  It is not the raw field `beta a`.

## 2. Reverse-KL chain rule and block composition

For fixed `S,y`, the parent conditional law on outside completions is

\[
p^V_{S,y}(z)=\frac{e^{\beta H_{S,y}(z)}}{K_{\beta,S}(y)}.
\]

The uniform outside mean of `H` is zero.  Consequently

\[
\boxed{
F^V_S(y)=D\!\left(U_{V\setminus S}\middle\|p^V_{S,y}\right).
}
\tag{R31.7}
\]

In particular, (R31.5) is the reverse-KL chain rule, because

\[
D\!\left(U_{\{i\}}\middle\|\nu_\beta(i\mid y)\right)
=\log\cosh b^V_{S,i}(y).
\tag{R31.8}
\]

For any ordering `i_1,...,i_k` of the outside set, let
`S_t=S union {i_1,...,i_t}`.  Iteration gives the order-independent identity

\[
\boxed{
F^V_S(y)=
\mathbb E_{Z\sim U_{V\setminus S}}
\sum_{t=0}^{k-1}
\log\cosh b^V_{S_t,i_{t+1}}(y,Z_{i_1},\ldots,Z_{i_t}).
}
\tag{R31.9}
\]

More generally, for `S subset R subset V`, let
`p^V_{R\setminus S\mid y}` be the parent marginal of the intermediate
block, with `V\setminus R` integrated out.  KL chain composition gives

\[
\boxed{
F^V_S(y)=
D\!\left(U_{R\setminus S}\middle\|p^V_{R\setminus S\mid y}\right)
+\mathbb E_{z\sim U_{R\setminus S}}F^V_R(y,z).
}
\tag{R31.10}
\]

This is associative and every displayed block charge is nonnegative.
There is a complementary internal-Gibbs formula.  If `F^R_S` is the
harmonic cost in the induced model on `R` and `mu^R(.|y)` is its Gibbs
extension kernel, then

\[
F^V_S=F^R_S+\log\mathbb E_{\mu^R(.\mid y)}e^{F^V_R}
=F^R_S+\mathbb E_{p^V(.\mid y)}F^V_R
-D\!\left(p^V(.\mid y)\middle\|\mu^R(.\mid y)\right).
\tag{R31.11}
\]

The negative KL in the second form is one reason that summing independent
one-deletion estimates is not an exact endpoint composition theorem.

## 3. Exact soft-selector response and cavity screening

For `|S|=m`, abbreviate `F_S(d)=F^V_S(d[S])` and set

\[
\Phi(d)=-\log\mathbb E_{S\sim U_m}e^{-F_S(d)},
\qquad
\pi_d(S)=\frac{e^{-F_S(d)}}{\sum_Te^{-F_T(d)}}.
\tag{R31.12}
\]

The Gibbs variational formula is pointwise:

\[
\boxed{
\Phi(d)=\mathbb E_{S\sim\pi_d}F_S(d)
+D(\pi_d\|U_m).
}
\tag{R31.13}
\]

For any two full cuts `d,d'`, put
`delta_S=F_S(d'[S])-F_S(d[S])`.  Then

\[
\boxed{
\Phi(d')-\Phi(d)
=-\log\mathbb E_{S\sim\pi_d}e^{-\delta_S}.
}
\tag{R31.14}
\]

This is the exact finite-difference formula required by any Efron--Stein or
cavity argument.

For a vertex flip `d'=d^i`, `delta_S=0` when `i notin S`.  When `i in S`,
apply (R31.4) to `S_0=S\setminus{i}` and define the screening residual

\[
r^V_{S,i}(d)=b^V_{S_0,i}(d[S_0])-\beta a_{S_0,i}(d)
=\frac{F_S(d)-F_S(d^i)}2.
\tag{R31.15}
\]

Hence

\[
\boxed{
\Phi(d^i)-\Phi(d)
=-\log\mathbb E_{S\sim\pi_d}
\exp\{2\mathbf1_{\{i\in S\}}r^V_{S,i}(d)\}.
}
\tag{R31.16}
\]

The relevant endpoint coordinate response is therefore a posterior soft
average of the **screening residual** `b-beta a`.  Raw row-square control
sees `a`, not this residual, agreeing with the `A_6` higher-mode wall.

There is also exact freedom to remove selectorwise common levels.  For any
real constants `alpha_S`, let

\[
q_\alpha(S)=\frac{e^{-\alpha_S}}{\sum_Te^{-\alpha_T}},
\qquad
\Psi_\alpha(d)=-\log\mathbb E_{q_\alpha}
e^{-[F_S(d)-\alpha_S]}.
\]

Then

\[
\boxed{
\Phi(d)=\Psi_\alpha(d)
+\log\frac{\binom nm}{\sum_Se^{-\alpha_S}}.
}
\tag{R31.17}
\]

Thus `Phi` and `Psi_alpha` have identical fluctuations and produce the same
normalized endpoint likelihood.  Any composition estimate which charges
the uncentered levels `alpha_S` is intrinsically lossy.

## 4. The rigorous accumulation wall

From (10.899), with `P=f_beta nu_beta`,

\[
\operatorname{Ent}_{\nu_\beta}(f_\beta)
=-\mathbb E_P\Phi-\log\mathbb E_{\nu_\beta}e^{-\Phi}.
\]

Dropping the first term and applying Jensen twice gives the valid but crude
chain

\[
\boxed{
\operatorname{Ent}_{\nu_\beta}(f_\beta)
\le-\log\mathbb E_{\nu_\beta}e^{-\Phi}
\le\mathbb E_{\nu_\beta,S}F_S.
}
\tag{R31.18}
\]

By (R31.7), for each `S` the final quantity is exactly

\[
\mathbb E_{Y\sim(\nu_\beta)_S}F_S(Y)
=D\!\left((\nu_\beta)_S\otimes U_{S^c}\middle\|\nu_\beta\right).
\tag{R31.19}
\]

It pays the full KL cost of replacing the parent outside conditional by
uniform independent spins.  Equations (R31.9)--(R31.10) show that it is a
sum of nonnegative reverse-KL charges.  Endpoint entropy, in contrast,
depends only on the centered fluctuation in (R31.17).

This is not merely a formal warning.  For the exact order-six minimizer
`A_6`, at `beta=1/2` and fixed density `m=3`, exact enumeration gives

\[
\boxed{
\operatorname{Ent}(f)=0.001953522524668\ldots,
\quad
-\log\mathbb E_\nu U=2.228571108927248\ldots,
\quad
\mathbb E_{\nu,S}F_S=2.774144107250946\ldots.
}
\tag{R31.20}
\]

Thus the two uncentered relaxations lose more than three orders of magnitude
already in an exact minimizer.  This is a finite mechanism wall, not an
asymptotic counterexample to (10.800).

The scaling audit reaches the same conclusion.  At
`beta=Theta(n^(-1/2+c))`, exact minimality gives
`||A||_op=O(n^(3/4))`, so the proved one-deletion charge is only
`O(beta ||A||_op)=O(n^(1/4+c))`.  Charging it at `k=Theta(n)` separate
steps would cost `O(n^(5/4+c))`, versus the required
`O(n^(1/2-2c))`.  The exact block identity does not remove this charge; it
shows that a proof has to center/cancel it.

## 5. Doob, Efron--Stein, and the sharp remaining criterion

Parameterize an oriented cut by its `n` independent sign bits (one
orientation bit and `n-1` vertex bits).  For

\[
\nu_{\beta,t}(d)\propto e^{-t\Phi(d)}\nu_\beta(d)
\]

and any bit-revealing filtration `F_j`, put
`M_{t,j}=E_{nu_(beta,t)}[Phi|F_j]`.  Doob orthogonality and (10.900) give
the exact identity

\[
\boxed{
\operatorname{Ent}_{\nu_\beta}(f_\beta)
=\int_0^1t\sum_{j=1}^{n}
\mathbb E_{\nu_{\beta,t}}(M_{t,j}-M_{t,j-1})^2\,dt.
}
\tag{R31.21}
\]

This is the sharp centered composition criterion.  A convergence-relevant
advance must prove, for target-specific exact minimizers at fixed density,

\[
\boxed{
\int_0^1t\sum_{j=1}^{n}
\mathbb E_{\nu_{\beta,t}}(M_{t,j}-M_{t,j-1})^2\,dt
=O(n^{1/2-2c}),
}
\tag{R31.22}
\]

with `Phi` evaluated through (R31.9), or equivalently through the posterior
screening responses (R31.16).  Unlike a bound on `E F`, (R31.22) discards
all selectorwise common accumulation.  It is an exact remaining lemma, not
a proof of it.

Ordinary Efron--Stein does not supply (R31.22).  First,
`nu_(beta,t)` is not a product law; a Poincare/tensorization comparison is a
new minimizer-specific input, and generic parent-Gibbs comparisons are
already blocked by (10.813).  Second, even if one granted a constant
comparison, the available pointwise response is far off scale.  If `i in S`,
flipping the child vertex changes each outside energy by at most `4k`, so

\[
|F_S(d^i)-F_S(d)|\le4\beta k,
\qquad |\Phi(d^i)-\Phi(d)|\le4\beta k.
\tag{R31.23}
\]

At fixed density the resulting sum of `n` squared bounds is
`O(n beta^2 k^2)=O(n^(2+2c))`, whereas (R31.22) requires
`O(n^(1/2-2c))`.  Exact minimality does not currently control the
renormalized screening residual in (R31.16), so replacing the raw bound by a
target-scale response is precisely the missing mathematical input.

A stronger but simpler sufficient statement, useful for falsification and
construction, is centered harmonic flatness: find constants `alpha_S` such
that

\[
\operatorname{osc}_{d}\Psi_\alpha(d)=O(n^{1/2-2c}).
\tag{R31.24}
\]

Then `Ent(f)<=osc Phi=osc Psi_alpha` proves the parent half of (10.800).
Equation (R31.9) makes (R31.24) a concrete assertion about the fluctuation,
not the total size, of accumulated effective-cavity reverse KL.  No such
uniform exact-minimizer theorem is presently known.

## 6. Frontier

Sequential integration does not prove the fixed-density parent estimate.
It does identify the correct object more sharply: the harmonic cost is an
accumulated reverse conditional KL, but only its centered selector-softened
screening response survives endpoint normalization.  First-moment,
uncentered entropy-chain, and naive block summation should not remain active.
The viable continuation is a minimizer-specific proof of either the exact
Doob energy (R31.22), preferably via the posterior screening formula
(R31.16), or the stronger centered flatness statement (R31.24).  The finite
`A_6` audit is a mechanism wall only; no convergence proof or asymptotic
falsifier is obtained.
