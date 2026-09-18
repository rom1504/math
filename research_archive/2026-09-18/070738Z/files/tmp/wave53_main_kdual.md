# Wave 53 main-agent scratch: an exact dual lower bound for the K profile

This is an independent derivation for later audit against the K-profile agent.

For `t=sqrt(H)`, infimal convolution duality gives the exact formula

```math
K_{1,2}(\beta,t)
=\sup\{\langle\beta,z\rangle:\lVert z\rVert_\infty\le1,
\ \lVert z\rVert_2\le t\}.
```

Indeed the conjugates of `||.||_1` and `t||.||_2` are the indicators of
the two displayed dual balls, and the conjugate of their infimal convolution
is their sum. For

```math
E=\lVert\beta\rVert_2^2,\qquad M=\lVert\beta\rVert_\infty,
```

take `z=alpha beta`, where

```math
\alpha=\min\{M^{-1},\sqrt{H/E}\}
```

(with the zero-vector case understood separately). This proves the exact
elementary lower bound

```math
\boxed{
\mathcal K_H(\beta)\ge
\min\left\{\frac{E}{M},\sqrt{HE}\right\}.}
\tag{M53.1}
```

Consequently, if `R=r+b_H`, then

```math
\boxed{
E\ge\max\left\{2MR,\frac{4R^2}{H}\right\}
\quad\Longrightarrow\quad
R\le\frac12\mathcal K_H(\beta).}
\tag{M53.2}
```

In the signing application `beta=2 sigma A[T,S]y`, hence

```math
E=4\lVert A[T,S]y\rVert_2^2,
\qquad M\le2m.
```

Thus the completely scalar sufficient condition

```math
\lVert A[T,S]y\rVert_2^2
\ge\max\left\{m(r+b_H),\frac{(r+b_H)^2}{H}\right\}
\tag{M53.3}
```

implies the dependence-safe tail (10.1276). At
`H=n^(3/4-c)` and `r=omega_n n^(3/2-c)`, the two requirements have scales
`omega_n n^(5/2-c)` and `omega_n^2 n^(9/4-c)`. This is feasible under the
known `O(n^(5/2))` cross-energy ceiling only for a sufficiently slow
`omega_n=o(n^c)` in the worst coordinate-concentrated case. It is stronger
than the necessary high-cross condition (10.1278), but it replaces the full
Holmstedt profile by two checkable scalars.

For a maximal-selector ground one can additionally combine (10.1271) with
`beta_j=2 sum_{u in S}a_{ju}y_u`: for every `i`,
`|beta_j|/2<=r_i+1`, so `M<=2(1+q_*/m)=O(sqrt(n))`. This observation does
**not** directly apply to the far-negative local states in (10.1277), which
need not be child grounds; it must not be used as an abundance proof.

The entropy-matched specialization is cleaner. Write

```math
X=\lVert A[T,S]y\rVert_2^2,
\qquad J=\lVert A[T,S]y\rVert_\infty,
\qquad R=r+b_H.
```

Since `E=4X` and `M=2J`, (M53.2) shows

```math
\boxed{
X\ge R^2/H,\qquad J\le R/H
\quad\Longrightarrow\quad
R\le\mathcal K_H(\beta)/2.}
\tag{M53.4}
```

Thus the natural missing coefficient regularity is not merely high cross
energy: it is high cross energy without a coordinate larger than its
entropy-matched scale `R/H`. For `R` of order
`omega_n n^(3/2-c)`, this coordinate scale is
`omega_n n^(3/4)`, while `X` has scale
`omega_n^2 n^(9/4-c)`.

There is also a useful scale audit. For uniform local `y`, Hanson--Wright on
`C=A[S,T]A[T,S]`, with `||C||op=O(n^(3/2))` and
`||C||F=O(n^(7/4))`, gives the upper bound `exp{-Omega(H)}` at deviations
`n^(9/4-c)`. In contrast, using only the trivial `J<=m` in (M53.3) demands
cross energy `n^(5/2-c)`, whose generic upper-tail exponent is at least
order `n^(1-c)`, much larger than `H=n^(3/4-c)`. Hence an abundance proof
based on (M53.2) must obtain genuine coordinate delocalization (or use the
full truncated profile); the trivial coordinate cap is entropy-mismatched.
