# Wave 51: direct and higher-moment completion tails

## Status

- **Verified:** the involution `w -> -w` gives an exact two-point formula for
  the completion CDF.  Combining it with Paley--Zygmund for the linear part
  and Cantelli for the outside quadratic part gives an explicit
  linear-dominance lower-tail bound.
- **Verified:** Bonami hypercontractivity gives a universal constant lower
  tail in the narrow layer `r < sqrt(V)/18`.
- **Verified:** a cubic endpoint polynomial gives an exact third-moment
  strengthening of the Wave 50 support/variance bound, and its third moment
  has a closed matrix formula.
- **Falsified as an independent far-tail mechanism:** positivity of the cubic
  numerator always implies positivity of the old quadratic numerator.  Thus
  the cubic bound can improve constants but can never enter a state on which
  the Wave 50 reverse branch vanishes.
- **Scoped obstruction:** the Bonami and constant-probability
  linear-dominance mechanisms require `r=O(sqrt(V))` or `r=O(s)`.  The existing
  exact-minimizer operator estimate gives both scales `O(n^(5/4))`.  For every
  `c<1/4` this lies inside the `O(n^(3/2-c))` local-margin tolerance, so saved
  mass of these states already gives the local restriction recurrence.  These
  mechanisms do not bypass the inverse/local-margin wall.
- **Numerical:** exhaustive checks on all sixteen stored matrix/selector
  pairs verify every identity and lower bound.  The cubic bound improves the
  prior envelope inside its support.  The pairing bound also detects states
  outside the two-moment support in three high-density finite cases, but these
  have one outside spin and do not defeat the asymptotic `n^(5/4)` wall.
- **Open:** a genuinely far-negative completion theorem would need a
  moderate-deviation lower bound beyond one conditional standard deviation,
  using coefficient regularity or a signed spectral property specific to
  exact minimizers.  Generic degree-two moment information does not provide
  it.

## 1. Notation and the exact pairing identity

Fix `S`, an oriented local state `(sigma,y)`, and `T=S^c`, `k=|T|`.  Write

```math
Z(w)=L(w)+Q_T(w),
\qquad
L(w)=\beta^{\mathsf T}w,
\qquad
Q_T(w)=w^{\mathsf T}Bw,
```

where

```math
\beta=2\sigma A[T,S]y,
\qquad B=\sigma A[T].
```

Both components have mean zero.  Put

```math
s^2=\mathbb E L^2=\lVert\beta\rVert_2^2,
\qquad
W=\mathbb E Q_T^2=2k(k-1),
\qquad
V=s^2+W.
```

Let `[w]` be uniform on the projective pairs `{w,-w}`.  Since `Q_T` is even
and `L` is odd, for every real threshold `u` there is the **Verified exact
identity**

```math
\boxed{
\Pr_w\{Z(w)\le u\}
=\mathbb E_{[w]}\frac12\left[
\mathbf 1_{\{Q_T(w)-|L(w)|\le u\}}
+\mathbf 1_{\{Q_T(w)+|L(w)|\le u\}}
\right].}
\tag{R51N.1}
```

In particular,

```math
\Pr\{Z\le-r\}\ge
\frac12\Pr\{|L|-Q_T\ge r\}.
\tag{R51N.2}
```

This is an exact use of the linear-plus-quadratic structure, rather than a
generic upper-tail estimate.

## 2. A verified linear-dominance lower bound

For `0<theta<1`, the fourth moment of the linear Rademacher sum obeys

```math
\mathbb E L^4
=3s^4-2\sum_i\beta_i^4\le3s^4.
```

Paley--Zygmund applied to `L^2` therefore gives

```math
\Pr\{|L|\ge\theta s\}\ge\frac{(1-\theta^2)^2}{3}.
```

If `v=theta*s-r>0`, one-sided Cantelli gives

```math
\Pr\{Q_T>v\}\le\frac{W}{W+v^2}.
```

The event `|L|>=theta*s` and `Q_T<=v` is contained in
`{|L|-Q_T>=r}`.  A union bound in (R51N.2) proves the **Verified theorem**

```math
\boxed{
\Pr\{Z\le-r\}\ge
\frac12\left[
\frac{(1-\theta^2)^2}{3}
-\frac{W}{W+(\theta s-r)^2}
\right]_+,\qquad \theta s>r.}
\tag{R51N.3}
```

Thus a fixed positive bound results, for example, when `r` is a fixed
subfraction of `s` and `W/s^2` is sufficiently small.  This can be positive
even when the Wave 50 support/variance numerator `V-r(q-e)` is nonpositive.

The scope is nevertheless limited at the convergence scale.  The established
exact-minimizer spectral estimate `||A||_op=O(n^(3/4))` gives

```math
s=2\lVert A[T,S]y\rVert_2
\le2\lVert A\rVert_{\rm op}\lVert y\rVert_2
=O(n^{5/4}).
\tag{R51N.4}
```

Consequently the hypothesis `r<theta*s` in (R51N.3) forces
`r=O(n^(5/4))`.  Since `5/4 <= 3/2-c` for `c<=1/4`, saved population of
these states already falls inside the local-margin restriction tolerance.
The theorem is a real conditional CDF bound, but not an independent route
around the local recurrence.

## 3. The universal Bonami boundary layer

For every centered degree-at-most-two Rademacher polynomial, Bonami gives

```math
\lVert Z\rVert_4\le3\lVert Z\rVert_2.
```

Interpolation between `L^1`, `L^2`, and `L^4` yields

```math
\mathbb E|Z|
\ge\frac{V^{3/2}}{(\mathbb E Z^4)^{1/2}}
\ge\frac{\sqrt V}{9}.
```

Since `E Z=0`, `E Z_-=(1/2)E|Z|>=sqrt(V)/18`.  On the other hand,
Cauchy--Schwarz gives

```math
\mathbb E Z_-
\le r+\mathbb E[Z_-\mathbf1_{\{Z_->r\}}]
\le r+\sqrt V\Pr\{Z<-r\}^{1/2}.
```

Hence, for `V>0`, the **Verified boundary-layer estimate** is

```math
\boxed{
\Pr\{Z\le-r\}\ge\Pr\{Z<-r\}
\ge\left(\frac1{18}-\frac r{\sqrt V}\right)_+^2.}
\tag{R51N.5}
```

At a nonnegative threshold this also gives
`Pr{Z<=g/p_2}>=Pr{Z<0}>=1/324`.  More generally, any state with

```math
g\ge-\eta p_2\sqrt V,\qquad 0\le\eta<1/18,
```

has conditional favorable probability at least `(1/18-eta)^2` (using
`1/324` when `g>=0`).  This slightly enlarges the positive-margin
population.  But (R51N.4), together with `W=O(n^2)`, also gives
`sqrt(V)=O(n^(5/4))`; the same local-margin scope obstruction applies.

## 4. A cubic endpoint-polynomial bound

Retain Wave 50 notation

```math
e=\mathbb E(e+Z),\quad
a=q+e,\quad b=q-e,\quad
D=q^2-e^2-V,
```

so `Z in [-a,b]`.  For `0<=r<a`, consider

```math
P_3(z)=(z+a)(z+r)(z-b).
```

It is nonnegative on `[-a,-r]` and nonpositive on `[-r,b]`.  Put

```math
M_3(a,b,r)=\max_{r\le u\le a}(a-u)(u-r)(u+b).
\tag{R51N.6}
```

The unique interior maximizer is

```math
u_*=\frac{a-b+r+
\sqrt{(a-b+r)^2+3\{b(a+r)-ar\}}}{3}.
```

Writing `mu_3=E Z^3`, direct expansion gives

```math
\mathbb E P_3(Z)
=\mu_3+(a-b+r)V-abr
=\mu_3+2eV-rD.
```

Therefore the **Verified cubic lower bound** is

```math
\boxed{
\Pr\{Z\le-r\}\ge
\frac{[\mu_3+2eV-rD]_+}{M_3(a,b,r)}.}
\tag{R51N.7}
```

A cruder denominator free of the maximization is
`M_3 <= (a-r)^2(a+b)/4`.

The third moment is itself explicit.  Walsh orthogonality and triangle
counting give

```math
\boxed{
\mu_3
=6\beta^{\mathsf T}B\beta+8\operatorname{tr}(B^3).}
\tag{R51N.8}
```

Indeed `3 E(L^2 Q_T)=6 beta^T B beta`, `E(L Q_T^2)=0`, and
`E Q_T^3=8 tr(B^3)`.

### Exact domination by the Wave 50 support

Let `phi_2(z)=(z+r)(z-b)`, whose expectation is `V-rb`.  On the tail
`z<=-r`, `phi_2>=0` and `z+a<=a-r`; off the tail, `phi_2<=0` and
`z+a>=a-r`.  Pointwise,

```math
(z+a)\phi_2(z)\le(a-r)\phi_2(z).
```

Taking expectations proves the **Verified decisive comparison**

```math
\boxed{
\mu_3+2eV-rD
\le(a-r)\{V-rb\}.}
\tag{R51N.9}
```

Thus a positive cubic numerator necessarily has `V>r(q-e)`, exactly the
positivity condition for the Wave 50 reverse two-moment branch.  The cubic
polynomial can sharpen the numerical lower bound because its maximum is
smaller after endpoint localization, but it cannot unlock even one new
conditional state.  In particular, it cannot escape the `g>=-O(n)` support
wall derived from the Wave 50 numerator in the relevant principal-gap
regime.

## 5. Exhaustive finite audit

`tmp/negative_completion_r51_check.py` exhausts every selector, oriented
local projective state, and outside completion for `A6`, `A8`, `A9`, and one
deterministic exact order-ten minimizer.  It verifies:

1. the variance and third-moment formulas;
2. the endpoint-polynomial expectation and all lower bounds;
3. the exact pairing identity in 33,600 separate tests;
4. Bonami's fourth-moment ceiling (the largest observed ratio is
   `E Z^4/V^2 = 5.2397`, below 81); and
5. that the combined bound never exceeds the exact CDF.

Selected averages over negative-margin conditional states are:

| case | exact CDF | old `H2` | cubic `H3` | pairing `HLQ` | max combined |
|---|---:|---:|---:|---:|---:|
| `A9,m=8` | 0.261663 | 0.005167 | 0.010370 | 0.037903 | 0.040778 |
| `A10,m=8` | 0.250000 | 0.008724 | 0.015052 | 0.000834 | 0.018860 |
| `A10,m=9` | 0.317164 | 0.008807 | 0.007656 | 0.089944 | 0.093712 |

No tested state has positive `H3` with zero `H2`, as (R51N.9) requires.
The pairing bound is positive while `H2=0` on `39.15%` of negative-margin
states for `A9,m=8` and `47.51%` for `A10,m=9`; these cases have `k=1`, so
the outside quadratic noise vanishes.  They validate the distinct mechanism
in (R51N.3), but are not fixed-density asymptotic evidence.

Reproduction:

```bash
.venv/bin/python tmp/negative_completion_r51_check.py \
  | diff -u tmp/negative_completion_r51_check.out -
```

## 6. Research judgment

The cubic route should be retained only as a sharper diagnostic inside the
existing two-moment support and retired as a proposed far-tail mechanism.
The Bonami boundary layer and the constant-probability form of the pairing
route are also trapped inside the already sufficient `O(n^(5/4))`
local-margin scale.

A noncircular direct completion attack must therefore control genuinely
moderate deviations `r >> sqrt(V)` (or `r >> s`) with probability at least
`exp{-O(n^(3/4-c))}`.  Equation (R51N.1) identifies two possible exact
inputs: a regular-coefficient lower tail for `L`, or a sufficiently negative
spectral/chaos tail for `Q_T`.  Either needs a minimizer-specific profile
theorem; moments through order three, support, and generic degree-two
hypercontractivity do not supply it.
