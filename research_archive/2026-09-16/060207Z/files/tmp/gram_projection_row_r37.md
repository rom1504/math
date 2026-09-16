# Wave 37 main attack: Gram projection removes selector representativeness

## Status

The identities and entropy/projection transfer theorem below are proved.  The
result does not prove convergence: it replaces the two-line row package
(10.1033) by one stronger but cleaner favorable-label Gram-energy estimate.
It is useful because the rare selector family no longer has to be
degree-three representative.

## 1. Full-column local energy

For a full word `z` and selector `S`, define

```math
K_S(z)=\lVert A[:,S]z_S\rVert_2^2
=z_S^T(A^2)[S]z_S.
```

If `p=m/n` and `p_2=(m)_2/(n)_2`, direct diagonal/off-diagonal expansion
gives

```math
\mathbb E_{U_m}K_S(z)
=p_2R_2(z)+(p-p_2)n(n-1).
\tag{G37.1}
```

For a singleton anchor `v`, put

```math
\alpha_1=\frac{m-1}{n-1},
\qquad
\alpha_2=\frac{(m-1)_2}{(n-1)_2},
\qquad
T_v=z_v(A^2z)_v-(n-1).
```

Then exactly

```math
\mathbb E[K_S(z)\mid v\in S]
=m(n-1)+\alpha_2\{R_2(z)-n(n-1)\}
+2(\alpha_1-\alpha_2)T_v.
\tag{G37.2}
```

Unlike `L_S(z)=||A[S]z_S||^2`, `K_S` sees the external rows.  The
conference-core/hub obstruction from Wave 36 is therefore visible already
inside a hub-avoiding selector: its positive hub rows contribute to
`K_S` even though the hubs are not selected.

## 2. Projection transfer from an arbitrary rare anchored law

Write the signed columns

```math
g_i=z_iA[:,i],
\qquad
V=\sum_i g_i=Az,
\qquad
R=\lVert V\rVert_2^2=R_2(z).
```

Assume first `R>0`, set `u=V/sqrt(R)`, and put
`a_i=<u,g_i>`.  Then

```math
\sum_i a_i=\sqrt R,
\qquad
\sum_i a_i^2
=\frac{z^TA^4z}{R}
\le\lVert A\rVert_{op}^2,
\qquad
|a_v|\le\sqrt{n-1}.
\tag{G37.3}
```

For `S` containing `v`, define `Z_S=sum_(i in S)g_i`.  Its scalar
projection has anchored-slice mean

```math
\mathbb E_{U_v}\langle u,Z_S\rangle
=\alpha_1\sqrt R+(1-\alpha_1)a_v.
\tag{G37.4}
```

Independent Bernoulli sampling on the other `n-1` vertices, global
Hoeffding with variance proxy `sum_(i!=v)a_i^2`, and conditioning on the
modal sample size `m-1` show that this scalar is subgaussian with proxy
`O(||A||_(op)^2)` and conditioning defect `O(log n)`.  Hence entropy
duality gives, for every law `P` on the anchored slice and

```math
H_P=D(P||U_v)+O(\log n),
```

the lower bound

```math
\mathbb E_P\langle u,Z_S\rangle
\ge\alpha_1\sqrt R-(1-\alpha_1)\sqrt{n-1}
-C\lVert A\rVert_{op}\sqrt{H_P}.
\tag{G37.5}
```

On the other hand,

```math
\mathbb E_P\langle u,Z_S\rangle
\le\mathbb E_P\lVert Z_S\rVert_2
\le\sqrt{\mathbb E_PK_S(z)}.
\tag{G37.6}
```

At fixed selector density, (G37.5)--(G37.6), including the trivial `R=0`
case, prove the deterministic transfer theorem

```math
\boxed{
R_2(z)
\le C_\rho\left\{
\mathbb E_PK_S(z)
+\lVert A\rVert_{op}^2[D(P||U_v)+\log n]
+n
\right\}.
}
\tag{G37.7}
```

This theorem needs no comparison of `P` with the uniform slice on pair or
triple inclusion probabilities.  Entropy is used only to transport one
scalar projection aligned with the unknown full row vector `Az`.

## 3. Agreement transfer and the exact remaining lemma

Let `y^S` be anchored favorable labels, let `z` be their coordinatewise
plurality, and put `P=U(G)`.  If
`e_S=d_H(z_S,y^S)`, principal compression gives pointwise

```math
|\sqrt{K_S(z)}-\sqrt{K_S(y^S)}|
\le2\lVert A\rVert_{op}\sqrt{e_S}.
```

Minkowski and the exact decoder (10.1002) give

```math
\sqrt{\mathbb E_PK_S(z)}
\le\sqrt{\mathbb E_PK_S(y^S)}
+2\lVert A\rVert_{op}\sqrt{\mathcal C_{v_*}(\mathbf y)}.
\tag{G37.8}
```

Combining (G37.7)--(G37.8) yields

```math
\boxed{
R_2(z)\le C_\rho\left\{
\mathbb E_PK_S(y^S)
+\lVert A\rVert_{op}^2
  [\mathcal C_{v_*}(\mathbf y)+\log\beta^{-1}+\log n]
+n
\right\}.
}
\tag{G37.9}
```

At the Wave 34 scales, exact minimality gives
`||A||_(op)^2=O(n^(3/2))`, the conflict contribution is below the target,
and

```math
\lVert A\rVert_{op}^2\log\beta^{-1}
=O(n^{9/4-c'}).
```

Consequently the single new sufficient estimate

```math
\boxed{
\mathbb E_{S\sim U(G)}
\lVert A[:,S]y^S\rVert_2^2
=O(n^{9/4-c'})
}
\tag{G37.10}
```

for some affordable anchored family with project-scale conflict proves the
weaker arbitrary-cut row clause directly.  This replaces both favorable
local-row regularity and one-sided degree-three representativeness in
(10.1033).  To reach the stronger center cap `2n(n-1)`, every term in
(G37.9), including the entropy term, would need the corresponding stronger
scale, so (G37.10) is presently an arbitrary-cut implementation only.

## 4. Why the result is not yet a proof

Decomposing the Gram energy gives

```math
K_S(y)=\lVert A[S]y\rVert_2^2
+\lVert A[S^c,S]y\rVert_2^2.
\tag{G37.11}
```

Generic child-ground stability controls only the first term, and only at
`O(n^(5/2))`; it says nothing power-saving about the external term.  Parent
fiber width also does not automatically control (G37.11): the external term
is the squared Fourier level-one norm of the completion polynomial, while a
large completion width is a range lower bound.  Wave 37's separate paired
replacement attack should be used to decide whether exact minimality can
select favorable labels satisfying (G37.10).

The exact checker `tmp/gram_projection_row_r37_check.py` verifies
(G37.1)--(G37.4), the Hamming inequality, and the completion-average identity

```math
\mathbb E_{x_{S^c}}R_2(y,x)
=K_S(y)+(n-m)(n-1).
```

