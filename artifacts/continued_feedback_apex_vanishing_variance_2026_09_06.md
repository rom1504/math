# Actual bounded-op signings with vanishing cubic feedback variance

Date: 2026-09-06. Status: exact derivation and integer replay completed;
submitted for independent audit. This falsifies a uniform pointwise
positive lower bound for the nonlinear feedback variance. It is not a
family of purported minimizers.

## 1. Construction and operator bound

Let `m>=2` be a power of two, let `H_m` be the symmetric Sylvester
Hadamard, and let `d=diag(H_m)`. Then `H_m²=mI`, its entries are signs,
and `sum d=0`. Set

```math
R=\begin{pmatrix}1&-1\\-1&1\end{pmatrix},\qquad
A=\begin{pmatrix}0&\mathbf1^{\mathsf T}\\
\mathbf1&R\otimes H_m-I_2\otimes\operatorname{diag}(d)
\end{pmatrix},\qquad B=\frac{A}{\sqrt{2m}}.
```

This is an actual symmetric hollow signing of order `2m+1`.
Before deleting the bulk diagonal, the normalized twin block acts only
on the antisymmetric pair subspace, with operator norm `sqrt2`. The
apex star acts only on the apex and the symmetric all-ones bulk vector,
with operator norm one. These subspaces are orthogonal. The hollowing
correction has operator norm `1/sqrt(2m)`, so

```math
\|B\|_{\mathrm{op}}\le\sqrt2+\frac1{\sqrt{2m}}.
```

## 2. Exact vanishing-variance formula

Put `Q=B²` and `T=B(Q^{circ3})B`. Let coordinate zero be the apex.
The exact claim is

```math
T_{00}=\frac{36m^2-40m+1}{8m^3}\longrightarrow0.             (1)
```

Let `V=A²` and index bulk coordinates by `(a,u)` with `a in {1,2}`,
`u in {1,...,m}`. Write `r_1=1,r_2=−1`. Squaring the displayed block
matrix gives

```math
V_{(a,u),(b,v)}
=2m r_a r_b\mathbf1_{u=v}
 -r_a r_b H_{uv}(d_u+d_v)
 +\mathbf1_{a=b,u=v}+1.                                     (2)
```

Because the apex row of `A` is one on the bulk,

```math
T_{00}=\frac1{(2m)^4}\sum_{a,b,u,v}V_{(a,u),(b,v)}^3.          (3)
```

For `u=v`, the two terms `a=b` have value `2m`, and the two terms
`a!=b` have value `3−2m`. Summing these cubes over `u` gives

```math
m\,[2(2m)^3+2(3-2m)^3]=72m^3-108m^2+54m.
```

For `u!=v`, put `z=H_uv(d_u+d_v)`. Summing over `a,b` gives
`2(1−z)^3+2(1+z)^3=4+12z²`. Since `z²=2+2d_ud_v` and `sum d=0`,
the total over distinct `u,v` is

```math
28m(m-1)+24\bigl[(\sum_ud_u)^2-m\bigr]=28m^2-52m.
```

Adding these contributions yields `2m(36m²−40m+1)` in (3), proving
(1). The exact replay is in
`computations/continued_feedback_variance_falsifier_2026_09_06.py`.

## 3. What it rules out

For any fixed `L>sqrt2`, all sufficiently large members meet
`||B||op<=L`, while their minimum nonlinear cubic feedback variance
is at most the quantity in (1). Thus no positive constant depending
only on `L` can uniformly lower-bound every `T_ii` in this class.

Only one vanishing-variance root is needed for that falsification. This
does not disprove a positive lower bound outside a negligible root set,
or the already proved lower bound on the average standard deviation.
It also does not imply low Boolean cap: the construction is an
adversarial signing family, not an extremal candidate.

The smooth zero-first covariance theorem therefore genuinely benefits
from avoiding a minimum-variance assumption. A threshold theorem can
still be sought with an exceptional-small-variance-set estimate instead
of a pointwise floor.
