# Audit of the Wave 28 child entropy--reward criterion

## Conclusion

The general criterion (R28.E10) is not independent of the unknown mean
restriction bound at the claimed scale.  Under its entropy budget, the
largest possible mean child-energy reward is itself
`O(n^(3/2-c))`.  Consequently (R28.E10) forces

```math
\mathbb E_{S\sim U_m}Y_A(S)=O(n^{3/2-c}),
```

in the one-sided upper-bound sense needed by (10.794).  Conversely, that mean
bound implies (R28.E10) by choosing the uniform child law.  The *general*
criterion and the mean restriction lemma are therefore equivalent at this
scale.

The hard-tail condition (R28.E13) remains a legitimate sufficient
specialization, but it is stronger than the mean lemma rather than an
independent weakening of it.

## 1. Oriented-projective normalization

For a fixed selector `S`, the oriented-projective alphabet is

```math
\mathcal D_S=\{(\sigma,[y]):
\sigma\in\{\pm1\},\ [y]\in\{\pm1\}^m/(y\sim-y)\},
\qquad |\mathcal D_S|=2^m.
```

Lift `(\sigma,[y])` by choosing either representative `y` uniformly.
If `\mu_S` is a law on `\mathcal D_S`, its lifted law is

```math
\widetilde\mu_S(\sigma,y)=\frac12\mu_S(\sigma,[y]).
```

The uniform oriented-projective law lifts to the uniform law on the
`2^(m+1)` pairs `(\sigma,y)`.  The likelihood ratio is unchanged on every
two-point fiber, so

```math
D(\widetilde\mu_S\Vert U_{\sigma,y})
=D(\mu_S\Vert U_{\mathcal D_S})=K_S.
\tag{A28.1}
```

Thus lifting costs neither a bit nor a factor two in entropy.  Under the
uniform reference, `sigma` and `y` are independent Rademachers and

```math
c_S(\sigma,y)=\sigma y^{\mathsf T}A[S]y
```

is centered and symmetric.  Its normalization is the same as in the ledger:
the quadratic form counts both matrix orientations, and

```math
\lVert A[S]\rVert_F^2=m(m-1),\qquad
\lVert A[S]\rVert_{\rm op}\le\lVert A\rVert_{\rm op}.
\tag{A28.2}
```

The extra orientation sign does not enlarge the absolute tail:
`|c_S|=|y^T A[S]y|`.  It only symmetrizes it.

## 2. Entropy-dual reward bound

The Rademacher Hanson--Wright mgf, equivalently its standard subgamma form,
gives universal constants `C_0,c_0>0` such that

```math
\log\mathbb E_{U_{\sigma,y}}e^{\lambda c_S}
\le C_0\lambda^2\lVert A[S]\rVert_F^2,
\qquad
0\le\lambda\le\frac{c_0}{\lVert A[S]\rVert_{\rm op}}.
\tag{A28.3}
```

Entropy duality and (A28.1) give, for every admissible `lambda`,

```math
\mathbb E_{\mu_S}c_S
\le\frac{K_S}{\lambda}
+C_0\lambda\lVert A[S]\rVert_F^2.
\tag{A28.4}
```

Optimizing in the quadratic range and otherwise taking the endpoint yields

```math
\boxed{
\mathbb E_{\mu_S}c_S
\le C\left(
\lVert A[S]\rVert_F\sqrt{K_S}
+\lVert A[S]\rVert_{\rm op}K_S
\right).
}
\tag{A28.5}
```

For an exact minimizer, `\lVert A\rVert_{\rm op}^2\le2q_n=O(n^(3/2))`.
Hence, uniformly in `S`,

```math
\mathbb E_{\mu_S}c_S
\le C\left(n\sqrt{K_S}+n^{3/4}K_S\right).
\tag{A28.6}
```

Writing `K=\mathbb E_S K_S`, Jensen gives the averaged bound

```math
\boxed{
\mathbb E_{S,d_S}c_S(d_S)
\le C\left(n\sqrt K+n^{3/4}K\right).
}
\tag{A28.7}
```

No absolute value is needed: only an upper bound on the positive reward is
used.  Applying the same argument to `-c_S` would also bound its negative
mean if desired.

## 3. Comparison at the target exponents

If `K=O(n^(3/4-c))`, then

```math
n\sqrt K=O(n^{11/8-c/2}),\qquad
n^{3/4}K=O(n^{3/2-c}).
\tag{A28.8}
```

For `0<c\le1/4`,

```math
\frac{11}{8}-\frac c2\le\frac32-c,
```

so (A28.7) becomes

```math
\mathbb E c_S=O(n^{3/2-c}).
\tag{A28.9}
```

Let

```math
L=\mathbb E_{S,d_S}
[Y_A(S)-(1-p_2)c_S(d_S)].
```

The first line of (R28.E10) is used as the one-sided bound
`L\le Cn^(3/2-c)`; an absolute `O` statement is stronger and also suffices.
Since `p_2` is fixed for the target pair and `0\le1-p_2\le1`,

```math
\mathbb E_SY_A(S)
=L+(1-p_2)\mathbb E c_S
\le O(n^{3/2-c}).
\tag{A28.10}
```

This is exactly the missing mean optimized-restriction estimate

```math
\mathbb E_SQ(A[S])
\le p^{3/2}q_n+O(n^{3/2-c}).
\tag{A28.11}
```

It immediately supplies a favorable selector.

Conversely, suppose (A28.11) holds.  Choose `\mu_S=U_{\mathcal D_S}` for
every selector.  Then `K_S=0` and the independent uniform orientation gives
`\mathbb E c_S=0`.  Therefore the first line of (R28.E10) is precisely the
mean bound (A28.10), while its entropy line holds with zero.  This proves the
claimed scale equivalence of the general criterion.

Finally, (R28.E13) constructs particular nonuniform upper-level laws which
satisfy (R28.E10).  The preceding implication shows that (R28.E13) itself
forces (A28.11).  The converse need not hold: a small mean envelope gives no
lower bound on every selector's high-energy level-set mass.  Thus the
hard-tail formulation is a stronger sufficient specialization, not an
equivalent characterization.
