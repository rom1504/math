# A continuous bounded function missed by all polynomial moments

2026-09-06. Elementary director counterexample for the richer-frame
approximation audit. This is NOT a counterexample to the proved first-marked
signing theorem, and is not a realized rich signing return. It shows why a
generic finite-degree Gaussian polynomial assumption does not license the
polynomial approximation of all bounded continuous functions of that output.

Let X be standard Gaussian and W=X^3. Define the bounded continuous function

```math
h(w)=\sin\left(\frac{\sqrt3}{2}|w|^{2/3}-\frac\pi6\right).
```

Then h is nonzero in L2 of the law of W, but

```math
\mathbb E[W^k h(W)]=0\quad\text{for EVERY integer }k\ge0.   (1)
```

Odd k vanish by parity. For k=2m, change variables `u=x^2` in the
Gaussian integral. Up to a positive real factor the integral is

```math
\operatorname{Im}\left[
 e^{-i\pi/6}\Gamma(3m+1/2)
 (1/2-i\sqrt3/2)^{-(3m+1/2)}\right].
```

The complex number in parentheses is `exp(-i*pi/3)`. Using the branch
with positive real part in the convergent gamma integral, the displayed
phase equals `-pi/6+(3m+1/2)pi/3=m*pi`; its imaginary part is zero.
This proves (1). Continuity and h(0)=-1/2 imply E h(W)^2>0.

Consequently, for every real polynomial p,

```math
\mathbb E[h(W)-p(W)]^2
=\mathbb E h(W)^2+\mathbb E p(W)^2
\ge\mathbb E h(W)^2>0.                              (2)
```

This is stronger than a warning about moment growth: it provides an exact
bounded continuous response which no polynomial catalog can approximate,
even for a SINGLE fixed law with all moments finite. Moreover the two
probability measures `(1+h)dmu` and `(1-h)dmu`, where mu is the law of W,
have identical moments and are different; nonnegativity and normalization
follow from |h|<=1 and (1) at k=0.

Appending any fixed independent Gaussian vector G to the retained query
list does not help: h(W) remains orthogonal to every joint polynomial in
(G,W), by independence. Retaining X itself WOULD help, since polynomials
in the original Gaussian X are dense. The example therefore concerns the
information actually retained, not an assertion that every polynomial
Gaussian construction loses bounded-response approximability.

The first-marked proof uses a separately proved exponential-moment property
of its literal coherent list, so this example does not challenge that
theorem. For a richer list, an exact bounded-smooth/discrete-Stein argument
could bypass polynomial density altogether. That is a possible alternative,
not something ruled out by (1). In particular this artifact must not be cited
as a falsehood proof for richer Boolean feedback closure.
