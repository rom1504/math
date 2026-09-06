# Arbitrary creation masks: the zero-strip slack obstruction

Date: 2026-09-05/06. This is a scoped final-campaign audit, not a proof of
uniform escape from all nonnegative masks.

Let `U` be the already constructed isometry from jointly even Gaussian
`L²` into the first Gaussian chaos. For `0≤H≤1`, write

```math
W=UH,\qquad v=E H^2=E W^2,\qquad
J(H)=E|W|(1-H).
```

The open question here is whether `J(H)≥.43` forces constants `κ,t0>0`,
uniform over all such creation masks, for which

```math
E[(1-H)1_{\{|W|≤t\}}]≥κt\qquad(0<t≤t0).             \tag{1}
```

The established central-mask theorem proves the needed statement in its
own class. This note does **not** extend that conclusion to arbitrary
masks.

## 1. Exact counterexample if the high-value condition is removed

Define an even Lipschitz function

```math
f(z)=1-\frac14\min\{(|z|-1)_+,1\}.
```

Then `3/4≤f≤1`, `f=1` on `[-1,1]`, and its weak derivative satisfies
`|f'|≤1/4`. Consequently

```math
v:=Ef(G)^2≥9/16,\qquad Ef'(G)^2≤1/16<v.
```

The scalar Gaussian creation fixed-point theorem (valid in Gaussian
`W^(1,2)`) supplies a unit first-chaos `V` with

```math
V=U[f(V)/\sqrt v].
```

Set `H=f(V)`. Then `W=UH=√v V` exactly and

```math
E[(1-H)1_{\{|W|≤t\}}]=0\qquad(0<t≤\sqrt v).
```

This is a valid creation mask, not an arbitrarily prescribed joint law.
It disproves an unconditional zero-strip slack assertion. It is not a
high-value counterexample: indeed

```math
J(H)≤\tfrac14 E|G|<\tfrac14<.43.
```

Smoothing the two corners away from `[-1,1]`, if desired, preserves a
smaller exact plateau and the strict derivative inequality.

## 2. A uniform finite-input criterion, with explicit constants

Suppose `H` is measurable with respect to a finite Gaussian coordinate
family `X`, and decompose the Gaussian linear form

```math
W=M(X)+\sigma Z,
```

where `Z` is standard normal independent of `X`. Assume

```math
\sigma^2≥\eta>0,\qquad J(H)≥j>0.                    \tag{2}
```

Since `v≤1`, Cauchy--Schwarz and `(1-H)^2≤1-H` give

```math
E(1-H)≥J(H)^2/v≥j^2.
```

Also `EM²≤1`. Set `R=√2/j`. Then

```math
E[(1-H)1_{\{|M|≤R\}}]≥j^2-P(|M|>R)≥j^2/2.
```

For `|M|≤R`, `0<t≤1`, and `η≤σ²≤1`, the conditional Gaussian density
throughout `[-t,t]` is at least

```math
d_\eta=(2\pi)^{-1/2}\exp[-(R+1)^2/(2\eta)].
```

Conditioning on `X` therefore proves the exact uniform estimate

```math
E[(1-H)1_{\{|W|≤t\}}]≥j^2 d_\eta\,t
\qquad(0<t≤1).                                      \tag{3}
```

Thus a lower bound on the genuinely independent output variance is a
sufficient missing hypothesis. Nondegeneracy of each finite instance
without a uniform variance bound does not give uniform constants.
Enlarging the conditioning family can only decrease this residual
variance, so it must be tracked for a specified input representation.

The argument also applies to any finite-dimensional Gaussian linear
input space, not only to a subset of the distinguished tree coordinates.
For the finite-coordinate masks used in the tree construction, taking an
ancestor-closed coordinate set makes the independent output remainder
literal: all output coordinates outside the set are independent of it.
Its variance is a sum of squared omitted Hermite coefficients. There is
no claim here that `J≥.43` bounds this sum below uniformly over growing
families.

### High value does not imply uniform intrinsic innovation

The following exact counterexample was supplied by the director and
independently checked here. It refutes a uniform lower bound in (2), even
when the conditioning space is the *minimal* Gaussian linear input
space of the mask. It does not refute (1).

Take the banked mask `H0=1{|V|≤α}`, `W0=UH0`, with
`J(H0)>.430658` and nondegenerate Gaussian pair `(V,W0)`. Put

```math
f(w)=\frac{1+\cos w}{2},\qquad
H_\epsilon=(1-\epsilon)H0+\epsilon f(W0),
\qquad 0<\epsilon<1.
```

This is jointly even and valued in `[0,1]`. The isometry and the
two-Lipschitz bound for `J` give

```math
J(H_\epsilon)≥J(H0)-2\epsilon>.43
\qquad(0<\epsilon<.000329).
```

By linearity,

```math
UH_\epsilon=(1-\epsilon)W0+\epsilon U[f(W0)].
```

Its squared residual norm outside `span(V,W0)` is consequently at most
`ε²||f||₂²≤ε²`.

The input space is genuinely two-dimensional. In the nonsingular
two-dimensional Gaussian plane, the mask has nonzero jumps on the two
lines `V=±α`. If it were measurable with respect to one linear form,
its distributional derivative in the orthogonal direction would vanish.
The jump measures force that direction to be orthogonal to `V`; on the
open regions away from the jumps, the nonconstant cosine term then
forces it to be orthogonal to `W0` as well. Nonsingularity makes this
impossible for a nonzero direction. Since the mask already depends only
on this plane, conditioning arguments exclude a smaller input space
using extraneous Gaussian directions.

The two input forms can have infinite distinguished-tree support. This
is allowed in the Gaussian variational closure; finite tree-coordinate
approximations follow from the established density theorem. The point
is that even a minimal Gaussian input space has vanishing innovation
along high-value masks. Therefore criterion (3), on its own, cannot
yield a uniform high-value theorem merely by minimizing the input span.
In fact this example retains uniform slack: `1-Hε≥(1-ε)(1-H0)`, while
the Gaussian pairs `(V,UHε)` converge in covariance to the nonsingular
pair `(V,W0)`. Their densities are uniformly positive on a fixed compact
rectangle outside `|V|≤α` and around output zero. Thus it deliberately
separates failure of the innovation criterion from failure of (1).

The same failure holds with **finite actual tree-coordinate inputs**.
First choose a fixed finite first-chaos form `V0` approximating the banked
`V` closely enough that `H0=1{|V0|≤α}` has `J(H0)>.43`; continuity of
thresholds in a nondegenerate Gaussian law and of `J` justifies this.
Write again `W0=UH0`. Choose increasing finite ancestor-closed coordinate
sets containing the support of `V0`, with projections `Pj`, and put
`Wj=PjW0`. Then `Wj→W0` in `L²`. Define

```math
H_{\epsilon,j}=(1-\epsilon)H0+\epsilon f(Wj).
```

For sufficiently large `j`, `(V0,Wj)` is nondegenerate. (The limiting
pair is nondegenerate: otherwise `H0` would be a scalar multiple of
`U^{-1}V0`, a bounded nonconstant finite polynomial, which is impossible.)
The preceding jump/gradient argument makes `span(V0,Wj)` the minimal
Gaussian input space, and the mask uses only the finite actual
coordinates in `Pj`. Since this span contains `Wj`, the residual norm of
`UH_(ε,j)` off either this span or the larger `Pj` is at most

```math
(1-\epsilon)\|W0-Wj\|_2+\epsilon\|Uf(Wj)\|_2
≤\|W0-Wj\|_2+\epsilon.
```

Choose `j` so that the first term is at most `ε`. The residual variance
is then at most `4ε²`, whereas
`J(H_(ε,j))≥J(H0)-2ε>.43` for sufficiently small `ε`. Thus the failure
of a uniform intrinsic innovation bound is not due to allowing infinite
Gaussian input forms.

## 3. Why a strip-filling feedback idea is not yet a counterexample

Let `H0` be a binary high-value mask and `W0=UH0`. For an even function
`q∈[0,1]`, equal to one on a central interval, the natural candidate is

```math
H=H0+(1-H0)q(W),\qquad W=UH.                        \tag{4}
```

Any actual solution would have an exact zero-slack strip. Put
`K=(1-H0)q(W)` and `Δ=UK`. The isometry and disjoint mask supports give
the useful exact identities

```math
W=W0+Δ,\qquad E[W0Δ]=E[H0K]=0,\qquad
EΔ²=EK²,\qquad EW²=EH0+EK².                        \tag{5}
```

Thus `Δ` and `W0` are independent Gaussian linear forms. Nevertheless,
if `H0=1{|V|≤α}`, the quantity

```math
E[VΔ]=E[(U^{-1}V)K]
```

has no sign determined by (5). Its change affects the certificate at
first order. The regression bound alone is

```math
|E[VΔ]|≤\sqrt{\operatorname{Var}(V\mid W0)\,EK²}.
```

Similarly, an actual local contraction proof for (4) must bound the
operator on all first-chaos directions `D` with quadratic form

```math
E[(1-H0)q'(W)^2D^2],                                \tag{6}
```

not merely the scalar average `E[(1-H0)q'(W)^2]`.
Directions correlated with the rare unmasked conditional tail can
increase (6). The proposed high-value construction has not passed these
two checks. No numerical covariance iteration or heuristic is being
reported as an existence or lower-value proof.

## 4. Precise remaining statement

The proved results establish (1) uniformly for the old central-mask
class, and establish (3) for arbitrary masks under the extra uniform
independent-variance condition (2). A high-value valid creation mask
violating (1), or a theorem deriving (1) from `J≥.43` alone, remains
unproved here. This is an actual remaining obstruction to transferring
the existing localized unmarked perturbation uniformly to the entire
arbitrary-mask variational class.
