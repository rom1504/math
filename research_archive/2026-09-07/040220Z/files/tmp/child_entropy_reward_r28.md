# Wave 28 root fallback: child entropy--reward removes the separate row budget

## Status

The general channel identities and the entropy-to-row theorem below are
proved.  They give an exact sufficient criterion at the retained-deficit
interface: a selector-indexed child law with entropy deficit
`O(n^(3/4-c))` and enough **mean child energy** automatically has the
required full-cut mutual information, effective loss, and row-square scales
after uniform outside extension.  The final audit below shows that this is
**not an independent proof mechanism**: the same entropy budget and
Hanson--Wright inequality bound the mean child-energy reward by
`O(n^(3/2-c))`, so the criterion is equivalent at this scale to the
already sufficient mean restriction lemma.

The hard-tail specialization is weaker than requiring every completion of a
near-ground state to be pointwise favorable: its energy threshold is the
actual restriction excess divided by `1-p_2`.  No minimizer-specific
reverse-tail theorem proving the entropy premise is known.  At a macroscopic
threshold, Hanson--Wright forces an `Omega(n^(3/4))` entropy deficit, so
generic quadratic concentration still misses every power saving.

## 1. Arbitrary child channel

Fix an exact order-`n` minimizer `A`, a target `m`, and
`S in binom([n],m)`.  The oriented projective child-cut space
`\mathcal D_S` has size `2^m`.  Choose an arbitrary law `\mu_S` on
this space, then extend the underlying relative spins independently and
uniformly outside `S`.  Let `D` be the resulting full oriented cut and
put

```math
K_S=D_{\rm KL}(\mu_S\Vert U_{\mathcal D_S})
=m\log2-H(\mu_S),
\qquad K=\mathbb E_{U_m}K_S.
\tag{R28.E1}
```

The outside extension contributes exactly `(n-m)log2` conditional
entropy, while the full oriented projective alphabet has size `2^n`.
Consequently

```math
\boxed{I(S;D)\le K.}
\tag{R28.E2}
```

For a child state `d_S=(sigma,[y])`, write
`c_S(d_S)=sigma y^T A[S]y`.  Every term involving an outside uniform
spin has zero mean, so

```math
\mathbb E[\langle A,D\rangle\mid S,d_S]=c_S(d_S).
\tag{R28.E3}
```

Using (10.806), or expanding (10.792), gives the exact conditional mean

```math
\boxed{
\mathbb E[\widehat\ell(S,D)\mid S,d_S]
=Y_A(S)-(1-p_2)c_S(d_S).
}
\tag{R28.E4}
```

Thus no pointwise control of the parent completion energy is needed by the
adaptive information inequality.

## 2. Entropy controls the row square

Let

```math
B_S=A[:,S]^{\mathsf T}A[:,S],
\qquad B_{S,0}=B_S-(n-1)I_m.
```

After lifting the projective child spin by an independent global sign,
data processing gives a law `\bar\mu_S` on `\{\pm1\}^m` with
`D(\bar\mu_S\Vert U)\le K_S`.  Uniform outside extension gives exactly

```math
\mathbb E[R_2(D)\mid S,Y]
=(n-m)(n-1)+Y^{\mathsf T}B_SY.
\tag{R28.E5}
```

The entropy variational formula combined with the Rademacher
Hanson--Wright mgf yields, for an absolute constant `C`,

```math
\mathbb E_{\bar\mu_S}Y^{\mathsf T}B_SY
\le \operatorname{tr}B_S
+C\left(
\lVert B_{S,0}\rVert_F\sqrt{K_S}
+\lVert B_{S,0}\rVert_{\rm op}K_S
\right).
\tag{R28.E6}
```

For completeness, apply
`E_\nu Z<=lambda^{-1}(D(\nu||U)+log E_U e^{lambda Z})` to
`Z=Y^T B_{S,0}Y`.  The quadratic-chaos mgf is
`O(lambda^2||B_{S,0}||_F^2)` while
`|lambda|<=c/||B_{S,0}||op`; optimizing `lambda` proves (R28.E6).

Here

```math
\begin{aligned}
\operatorname{tr}B_S&=m(n-1),\\
\lVert B_S\rVert_{\rm op}
&\le\lVert A\rVert_{\rm op}^2=O(n^{3/2}),\\
\lVert B_{S,0}\rVert_{\rm op}&=O(n^{3/2}),\\
\lVert B_{S,0}\rVert_F^2
&\le\operatorname{tr}(B_S^2)
\le\lVert B_S\rVert_{\rm op}\operatorname{tr}B_S
=O(n^{7/2}).
\end{aligned}
\tag{R28.E7}
```

Averaging (R28.E5)--(R28.E7) and using Jensen therefore proves

```math
\boxed{
\mathbb E R_2(D)
\le n(n-1)
+O\left(n^{7/4}\sqrt K+n^{3/2}K\right).
}
\tag{R28.E8}
```

In particular, for every fixed `0<c<1/4`,

```math
K=O(n^{3/4-c})
\quad\Longrightarrow\quad
\mathbb E R_2(D)=O(n^{9/4-c}).
\tag{R28.E9}
```

The first error has exponent `17/8-c/2<=9/4-c` exactly when
`c<=1/4`.

## 3. Exact sufficient entropy--reward criterion

Equations (R28.E2), (R28.E4), and (R28.E9), substituted in (10.794) with
`lambda asymp n^(-3/4)`, prove the following.

**Criterion.**  Fix `\rho\in[1/2,1)`, `0<c<1/4`, and uniform
constants. If for every sufficiently large `n` and every integer
`m\in[\rho n,n)` there are a target-specific exact minimizer and child
laws `\{\mu_S\}` such that

```math
\boxed{
\begin{aligned}
\mathbb E_{S,d_S}
\left[Y_A(S)-(1-p_2)c_S(d_S)\right]
&=O(n^{3/2-c}),\\
\mathbb E_SK_S&=O(n^{3/4-c}),
\end{aligned}
}
\tag{R28.E10}
```

then the power-saving principal-restriction edge, and hence convergence,
follows.  The row-square premise is automatic.

This is a nonlinear selector-local rate--distortion problem.  It does not
assume the unknown mean restriction lemma: the child energy reward is
retained explicitly.

## 4. Hard-tail specialization and its exact scale wall

For each selector define

```math
a_S(t)=\frac{Y_A(S)-t}{1-p_2},
\qquad
f_S(t)=U_{\mathcal D_S}
\{d_S:c_S(d_S)\ge a_S(t)\}.
\tag{R28.E11}
```

Take `t\ge0`. The set is nonempty: averaging uniform outside
completions of an exact child ground proves `Q(A[S])\le q_n`, and then
`p_2\le p^{3/2}` gives
`p_2Q(A[S])\le p^{3/2}q_n+t`.
Take `\mu_S` uniform on this upper level set.  Then

```math
K_S=-\log f_S(t),
\qquad
\mathbb E[\widehat\ell\mid S]\le t.
\tag{R28.E12}
```

Consequently the concrete open lemma

```math
\boxed{
\mathbb E_{S\sim U_m}[-\log f_S(t)]
=O(n^{3/4-c}),
\qquad 0\le t=O(n^{3/2-c}),
}
\tag{R28.E13}
```

is sufficient for convergence.  Its threshold depends on the nonlinear
envelope `Y_A(S)`, so it is not the unaligned scalar sign tail retired in
Wave 26.

It nevertheless has a sharp generic wall.  If
`a_S(t)>=gamma n^{3/2}` on a selector set of fixed positive mass, then
for each such selector

```math
\lVert A[S]\rVert_F^2=m(m-1)=O(n^2),
\qquad
\lVert A[S]\rVert_{\rm op}
\le\lVert A\rVert_{\rm op}=O(n^{3/4}).
```

Hanson--Wright gives

```math
f_S(t)
\le2\exp\{-\Omega(n^{3/4})\},
\tag{R28.E14}
```

and hence the left side of (R28.E13) is `Omega(n^(3/4))`, too large
for every fixed power saving.  At threshold
`a_S=O(n^(3/2-c))`, the operator-norm exponent becomes
`O(n^(3/4-c))`, so this obstruction is exactly scale-matched and no
longer rules out (R28.E13).

The missing positive input is therefore a minimizer-specific,
envelope-aligned high-energy reverse tail at the power-saving boundary.

## 5. Final circularity audit

The same entropy transport argument applied to the signed child energy
`c_S(\sigma,y)=\sigma y^{\mathsf T}A[S]y` gives

```math
\boxed{
\mathbb E_{S,d_S}c_S(d_S)
\le
O\left(n\sqrt K+n^{3/4}K\right).
}
\tag{R28.E15}
```

Indeed the uniform oriented projective law makes `c_S` centered, its
absolute tail is the Rademacher quadratic-form tail, and
`\lVert A[S]\rVert_F=O(n)`,
`\lVert A[S]\rVert_{\rm op}=O(n^{3/4})`. Apply entropy duality for each
selector and Jensen-average `\sqrt{K_S}`.

At `K=O(n^{3/4-c})`, with `0<c<1/4`, the right side of (R28.E15) is
`O(n^{3/2-c})`. Therefore (R28.E10) itself implies

```math
\mathbb E_{S\sim U_m}Y_A(S)=O(n^{3/2-c}),
\tag{R28.E16}
```

which is exactly the already sufficient mean restriction lemma (10.796).
Conversely, if (R28.E16) holds, take every child law uniform: then
`K=0`, the mean child energy is zero, and (R28.E10) holds.

Thus the general entropy--reward criterion is equivalent at the relevant
scale to the mean restriction lemma, while its hard-tail specialization is
at least as strong. Both are useful diagnostics for the `n^{3/4}` wall,
but should be retired as independent successors to compressed fixed-cut
alignment.
