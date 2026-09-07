# Independent reconstruction of the p=24/25, t=97/20 ternary upper

Status: **PASS**, mathematical reduction reconstructed and independent interval
cover replay completed. This proves an improved all-order original upper
bound through the existing direct-E realization. It does not prove convergence
or identify the actual optimum with the ternary construction.

```math
\limsup_n M_n/n^{3/2}
\le {t+p\log2-5151/6250\over2t\sqrt p}
<0.493608094,
\qquad p=24/25,\quad t=97/20.
```

The middle constant lies in the exact decimal-rational enclosure

```math
[0.493608093588748652719124428717,
 0.493608093588748652719124428718].
```

The earlier displayed floating approximation .49360809358874863 is not an
outward endpoint and is not used as an asserted inequality.

## 1. Rate-distortion and reproduction support, without a three-atom ansatz

For `nu_p=(1-p)delta_0+(p/2)(delta_{-1/sqrt(p)}+delta_{1/sqrt(p)})`, let
`J_lambda=inf_L[I(X;L)+lambda E Var(X|L)]`. Replacing a reproduction by its
posterior mean decreases squared error and cannot increase information.
Consequently reproduction values can be restricted to the compact source
hull. The Gibbs variational formula, together with minimizing the output
reference measure in relative entropy, gives

```math
J_\lambda=-\max_\pi\sum_x\nu_p(x)
 \log\int e^{-\lambda(x-y)^2}\,d\pi(y).
```

The objective is concave and continuous on reproduction measures; symmetry
can be imposed by averaging pi with its reflection. At a maximizer its KKT
witness is at most one everywhere and equals one on its support:

```math
q(y)=e^{-\lambda y^2}[a+b\cosh(2\lambda y/\sqrt p)],\qquad a,b>0.
```

Set r=a/b and z=2lambda y/sqrt(p). Positive critical points solve
`h(z)=p/(2lambda)`, where `h=sinh(z)/[z(r+cosh(z))]`. Direct differentiation
gives sign(h')=sign(r-R(z)), with

```math
R(z)={\sinh z\cosh z-z\over z\cosh z-\sinh z}.
```

The denominator is positive. Its coefficient at z^(2j+1) is
`2j/(2j+1)!`, while the numerator coefficient is `4^j/(2j+1)!`, j>=1.
The ratio `4^j/(2j)` is strictly increasing. Expanding the derivative of the
power-series quotient in pairs of unequal indices therefore proves that R
strictly increases from 2 to infinity. Hence h decreases, or first increases
then decreases. The derivative of q has the sign of h-p/(2lambda), so q has
at most one positive local maximum.

At the source endpoint y=1/sqrt(p), the positive-source Gaussian summand has
derivative zero while both other summands have strictly negative derivatives.
Outside the source hull all summands decrease. Thus this endpoint cannot be
a support maximum. Symmetry and the KKT equality prove that an optimal pi is
supported on `{0,+a0,-a0}`, with `0<=a0<=1/sqrt(p)`.

Every argument in this section applies to all 0<p<1 and lambda>0; there is
no hidden restriction to the previously certified p=31/32, t=4 phase.

## 2. Exact elimination of the central weight

Write pi=w delta_0+(1-w)(delta_a0+delta_-a0)/2 and z=a0 sqrt(p). Set

```math
A=e^{\lambda z^2/p}-1,\quad B=\cosh(2\lambda z/p)-1,
\quad v={(1-w)e^{-\lambda a0^2}\over w+(1-w)e^{-\lambda a0^2}}.
```

The change w->v is a continuous bijection of [0,1], reversing endpoints.
With K=w+(1-w)exp(-lambda a0²), one has `K(1+vA)=1`; the zero-source
kernel integral is K and each nonzero-source integral is
`exp(-lambda/p) K(1+vB)`. Therefore

```math
J_\lambda=\lambda-\max_{z,v\in[0,1]}G(A,B,v),
\qquad G=p\log(1+vB)-\log(1+vA).
```

For A,B>0, the sign of G_v is the sign of the strictly decreasing affine
numerator `pB-A-(1-p)ABv`. Thus its unique maximizer is

```math
v_*=\operatorname{clip}_{[0,1]}
 {pB-A\over(1-p)AB}
=\operatorname{clip}_{[0,1]}{p/A-1/B\over1-p}.
```

The endpoint cases A=0 or B=0 are elementary and are included in both
verifiers. G is decreasing in A and increasing in B for each fixed v,
hence its optimized value has the same monotonicities.

Using `g_t(v)=sup_{0<lambda<=t}[c_t(lambda)-lambda v]`, where
`c_t(lambda)=log(lambda(2t-lambda)/t²)/4`, the two suprema commute, giving

```math
E_t(\nu_p)=\sup_{0<\lambda\le t,\ 0\le z\le1}
 [c_t(\lambda)-\lambda+G(A,B,v_*)].                 (1)
```

## 3. Low precisions and the Gaussian comparison

For p>=1/3, the even moments p^(1-j) of nu_p are bounded by the Gaussian
moments (2j-1)!!, by induction; thus nu_p is unit subgaussian. Relative-entropy
duality at a posterior mean mu gives D(posterior||nu_p)>=mu²/2. Averaging
shows `J_lambda=lambda` for lambda<=1/2, the constant channel attaining the
lower bound. Our p=24/25 meets this hypothesis.

The global maximum of c_t(lambda)-lambda is g_t(1), at

```math
\lambda_*=t(1-\rho),\qquad
\rho={\sqrt{1+16t^2}-1\over4t}.
```

In fact lambda_*<1/4 for every t>0, so it lies inside the low-precision
region. A separate rational-only calculation gives

```math
g_t(1)\in[-0.824473344812113062358239490940,
           -0.824473344812113062358239490939],
```

strictly below the target -5151/6250=-0.82416. It also encloses lambda_*
between 0.243560975364768362878092138566 and the next 10^(-30) endpoint.

No claim is made that the new parameters are in the Gaussian phase: only
the upper bound on E in (1) is needed.

## 4. The envelope-gradient rectangle bound is valid at clipping transitions

On z>0, A,B>0 and the maximizing v is unique and continuous. The envelope
therefore has continuous derivatives given by differentiating at v_*,
including points where the optimizer reaches 0 or 1. There is no derivative
of v_* to pay. With `r=lambda z²/p` and `q=2lambda z/p`, these are

```math
F_\lambda={1\over4}\left({1\over\lambda}-{1\over2t-\lambda}\right)-1
 +{2v_*z\sinh q\over1+v_*B}
 -{v_*z^2e^r\over p(1+v_*A)},
```

```math
F_z={2v_*\lambda\sinh q\over1+v_*B}
 -{2v_*\lambda ze^r\over p(1+v_*A)}.
```

These match the proposed verifier exactly. On a closed rectangle with a>0,
an interval enclosure of the unclipped rational weight expression encloses
v_* after monotone clipping. Interval arithmetic then bounds both gradient
components. If their absolute bounds are L_lambda,L_z, the line segment from
the rectangle center proves

```math
F\le F(\text{center})+L_\lambda(u-l)/2+L_z(b-a)/2.
```

Rectangles touching z=0 do not use that formula. The alternate monotonic
bound is valid throughout: take A at (l,a), B at (u,b), and c_t(lambda)-lambda
at l. The latter function decreases on lambda>=1/2, because its derivative
is at most 1/(4lambda)-1<=-1/2.

## 5. Independent replay, not an imported PASS

The reconstruction program
`computations/flatify_adversary_2026_09_07_ternary_replay.py` imports no proposer code.
It differs in three substantive implementation choices:

- It evaluates G directly at the exact rational maximizing weight of the
  conservative A/B endpoints, avoiding the proposer's closed-form log branch.
- It bounds the optimizer using `(p/A-1/B)/(1-p)`, and factors common terms
  differently in both gradient formulas.
- It uses a different normalized-width split rule and 60 decimal digits.

The initial domain is exactly `[1/2,97/20] x [0,1]`. Rational midpoint splits
partition every rejected rectangle into its two exact children. All 21,538
accepted leaves have strict upper bounds below -5151/6250; there are 43,075
visited nodes and zero pending nodes. Of the leaves, 21,472 use the mean-value
bound. Exact area sums to 87/20, and the binary-tree node identity is checked.
The worst accepted bound is lower than the target by approximately 5.795e-9;
the actual comparison uses exact fractions. Run time was 83.24 seconds.

Every accepted leaf, its rational upper bound, and its method are saved in
`computations/results/flatify_adversary_2026_09_07_ternary_replay_leaves.json.gz`
(888,088 bytes). The summary JSON is in `computations/results/` with the same
base name as the script. The installed mpmath interval exp
and log source was directly inspected: lower endpoints use round_floor and
upper endpoints use round_ceiling. No floating approximation prunes a box.

The separate `computations/flatify_adversary_2026_09_07_rational_conversion.py` does not
use mpmath: logarithms are bracketed by an 80-term artanh series with an exact
geometric tail, and square roots by 40-digit integer-square-root brackets.
It checks the Gaussian comparison, lambda_*, positivity of the cap numerator,
and the outward cap enclosure stated at the start. Its exact output is saved
as `computations/results/flatify_adversary_2026_09_07_rational_conversion.json`.
The original `tmp/` scripts and outputs remain intact for normal archival
preservation; the canonical scripts now write directly to `computations/results/`.

## 6. Actual-sign realization at these exact parameters

The existing direct-E stopping proof fixes arbitrary t>0 and proves
`limsup_r B^r Phi_t(nu)<=E_t(nu)` using the Phi-G drift budget. Neither its
precision supersolution nor its stopping/continuity arguments require t=4
or Gaussian-phase equality. The fixed-depth row construction applies to any
fixed 0<p<=1: k=floor(pm), each normalized input has second moment exactly
one, and at fixed depth its finite alphabet depends only on p and that depth.

For every full spin and both polarities the exact defect identity is
`D_sigma=2(m²k-sigma x^T W_T x)`. The full-spin/fibre expectation estimate
has exponential rate `t gamma+p log2+B^r Phi_t(nu_p)`. Independent fibres
give `(E Z)^m`, not `E Z^m`. Therefore any
`gamma<-(p log2+E_t(nu_p))/t` has a favorable strict exponent after selecting
a sufficiently large but finite depth. At the resulting parent order N=mk,
the hollow energy cap is at most `(1-gamma)/(2sqrt(p))+o(1)`, including the
O(N) diagonal deletion. Let the strict margin vanish after the order limit.

This yields exactly `(t+p log2+E_upper)/(2t sqrt(p))`. Our rational check
verifies p log2+E_upper<0 and a positive numerator. H2/H12 terminal orders
remain multiplicatively dense with fixed depth, and principal restriction
fills all missing orders with factor 1+o(1). None of these steps depends on
the old p=31/32,t=4 values. No actual-optimal-child hypothesis, growing-depth
type theorem, or H=E identity is introduced by this parameter change.

No mathematical correction was found in the new certificate or its original
all-order consequence. The numeric improvement does not close the original
liminf-to-all-order gap.
