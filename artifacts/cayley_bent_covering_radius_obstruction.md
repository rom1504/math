# Cayley/bent covering-radius obstruction

Date: 2026-07-31. This is an agent-authored exact audit of a scale-dependent
coding/design mechanism. It treats arbitrary functions at every scale, not a
fixed tensor generator.

## 1. Cayley signing and exact Walsh diagonalization

Let `V=F_2^m`, `n=2^m`, and let

```math
g:V\longrightarrow\{+1,-1\}.
```

Define the translation-invariant signing

```math
(A_g)_{uv}=
\begin{cases}
g(u+v),&u\ne v,\\
0,&u=v.
\end{cases}                                           \tag{CB1}
```

This includes difference-set, bent-function, plateaued-function, and
first-order Reed--Muller concatenation constructions. No tensor or fixed
formula for `g` is assumed.

Write

```math
W_g(a)=\sum_{z\in V}g(z)(-1)^{a\mathbin\cdot z}       \tag{CB2}
```

for the Walsh transform, and put

```math
\chi_a(u)=(-1)^{a\mathbin\cdot u}.
```

The full Cayley matrix

```math
Q_g=A_g+g(0)I
```

is convolution by `g`. Therefore

```math
Q_g\chi_a=W_g(a)\chi_a,
\qquad
A_g\chi_a=(W_g(a)-g(0))\chi_a.                       \tag{CB3}
```

The crucial point is that every character `chi_a` is itself Boolean.

## 2. Exact cap theorem

For every Boolean `x`, the spectral estimate gives

```math
|H_{A_g}(x)|le {n\over2}
 \max_a|W_g(a)-g(0)|.                                 \tag{CB4}
```

Choose a character attaining the maximum in (CB4). Equation (CB3) gives
equality. Hence

```math
\boxed{
\operatorname{cap}(A_g)
={n\over2}\max_{a\in V}|W_g(a)-g(0)|.}               \tag{CB5}
```

This is the exact two-sided absolute quadratic normalization. There is no
one-sided maximum-cut substitution: the absolute value in (CB5) is attained
by a Boolean eigencharacter of the actual symmetric signing.

Parseval says

```math
\sum_{a\in V}W_g(a)^2=n^2.                            \tag{CB6}
```

Thus `max_a |W_g(a)|>=sqrt(n)`, and (CB5) yields the uniform theorem

```math
\boxed{
\operatorname{cap}(A_g)
\ge {n\over2}(\sqrt n-1).}                            \tag{CB7}
```

In particular, for every scale-dependent sequence `g_m`,

```math
\liminf_{m\to\infty}
{\operatorname{cap}(A_{g_m})\over n^{3/2}}
\ge\frac12.                                           \tag{CB8}
```

## 3. Exact covering-radius mapping

Write `g=(-1)^f` for a Boolean function `f:V->F_2`. The Hamming distance to
the affine word `a dot z+b` is

```math
d\bigl(f,a\mathbin\cdot z+b\bigr)
={n-(-1)^bW_g(a)\over2}.                              \tag{CB9}
```

Minimizing over the affine offset `b` gives

```math
d\bigl(f,\mathrm{RM}(1,m)\bigr)
={n-\max_a|W_g(a)|\over2}.                            \tag{CB10}
```

Therefore the covering radius of the first-order Reed--Muller code is

```math
\rho(\mathrm{RM}(1,m))
={n-\min_g\max_a|W_g(a)|\over2}.                      \tag{CB11}
```

Equations (CB5) and (CB10) differ only by the diagonal correction `g(0)`,
which changes the relevant maximum by at most one. Thus a theorem producing a
deep hole of `RM(1,m)` maps exactly to a Cayley signing with cap

```math
{n\over2}\bigl(n-2d(f,\mathrm{RM}(1,m))+O(1)\bigr),  \tag{CB12}
```

not to a smaller or one-sided graph quantity.

For even `m`, bent functions have

```math
|W_g(a)|=\sqrt n\quad\text{for every }a,              \tag{CB13}
```

so they attain the Parseval floor in (CB6). But (CB5) then gives normalized
cap `1/2+o(1)`, not a constant below `1/2`. For odd `m`, any plateaued or
concatenated replacement still obeys (CB7), and spectral non-flatness can
only increase the cap certificate.

## 4. Consequence for scale-dependent fusion

The obstruction is not the fixed-fiber persistence theorem in different
notation. One may choose `g_m` independently at every order, use a different
difference set at every scale, or apply any known nonlinear bent/plateaued
concatenation. Equation (CB7) survives because every resulting Cayley algebra
has Boolean extremal characters and Parseval forces one Walsh coefficient of
size at least `sqrt(n)`.

This gives a precise negative answer for the natural proposed mechanism:

> A first-order Reed--Muller covering-radius concatenation, mapped through an
> additive Cayley signing, cannot yield a structured family whose normalized
> cap is uniformly below `1/2`.

The native orders are powers of two, so a separate order-filling operation
would still be required. But such an operation cannot repair the landing
behavior on those native orders. If the true limiting constant is below
`1/2`, (CB8) gives a linear `b`-scale gap; if it equals `1/2`, proving landing
again requires the unknown sharp lower constant.

## 5. Surviving coding escape

A coding construction must break at least one exact hypothesis behind
(CB5):

1. use a non-translation association scheme whose extremal eigenvectors are
   not Boolean;
2. use a nonabelian fusion in which character dimensions grow and no Boolean
   one-dimensional channel carries the spectral maximum; or
3. prove a state-dependent amalgamation inequality that controls the Boolean
   cap without passing through the complete Walsh spectrum.

Merely replacing bent functions by a scale-dependent concatenated deep hole
does not qualify: (CB5)--(CB12) show that this is an exact reformulation whose
best possible constant is the existing conference value `1/2`.
