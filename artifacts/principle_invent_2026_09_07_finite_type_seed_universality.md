# Finite-type seed universality through a canonical dual

Date: 2026-09-07. Status: **Pending independent audit**. This extends the
symmetric-type pressure proof to heterogeneous, possibly asymmetric local
types. The signing may influence a finite partition, but its entire leading
pressure contribution is controlled by its bilinear Boolean norm.

## 1. Definitions and theorem

Fix a finite involutive alphabet E of size L, and a strictly positive kernel
K with `K(a,b)=K(b,a)=K(-a,-b)`. Put

```math
K_+(a,b)=K(a,b),\qquad K_-(a,b)=K(a,-b),\qquad
\kappa=\min K/\max K>0.
```

Both K_+ and K_- are symmetric. At each of m vertices choose an arbitrary
probability type `nu_i` on E, with all its positive masses at least a fixed
alpha>0, and `d nu_i(a)` integral for d=m-1. Zero masses may be handled by
restricting the local alphabet. Let the rows independently be uniform
permutations of these prescribed types and define

```math
Z_S=\mathbb E\prod_{i<j}K_{S_{ij}}(u_{ij},u_{ji}),\qquad
\beta(B)=\max_{x,y\in\{-1,1\}^m}|x^TBy|.
```

For simplex vectors r_i on the support of nu_i, define

```math
\Psi_S(r)=\sum_{i<j}\log(r_i^T K_{S_{ij}}r_j)
+d\sum_iD(\nu_i\Vert r_i).
\tag{1}
```

The following two assertions hold, with constants independent of m, the
types, and the signing:

```math
\log Z_S=\inf_r\Psi_S(r)+O_{L,\alpha,K}(m^{3/2}\log m),
\tag{2}
```

```math
\left|\inf_r\Psi_S(r)-\inf_r\Psi_T(r)\right|
\le C_{L,K}\,\beta(S-T).
\tag{3}
```

Consequently any two bounded-cap seed sequences, `Q(S),Q(T)=O(m^(3/2))`,
have the same leading m² kernel pressure, uniformly over all these
heterogeneous asymmetric types:

```math
|\log Z_S-\log Z_T|=O_{L,\alpha,K}(m^{3/2}\log m).
\tag{4}
```

The implication uses the elementary polarization bound `beta(S)<=4Q(S)`
for hollow symmetric S. This theorem compares actual annealed partitions.
It does not compare the actual energy maxima produced by two seeds.

## 2. Canonical dual and its positive optimizer

Let P_0 make all directed incidences independent, with row-i marginal nu_i.
Let T_types denote the prescribed exact-type event. Consider edge laws
`pi_ij` with the aggregate marginal constraints

```math
\sum_{j\ne i}\pi_{ij}^{(i)}(a)=d\nu_i(a).
```

The strictly concave entropy-energy optimization over these edge laws has
an interior feasible point (the product laws) and a strictly positive
unique maximizer. Lagrange multipliers h_i(a), with additive constants
irrelevant, give

```math
\pi_{ij}(a,b)=
\frac{\nu_i(a)\nu_j(b)e^{h_i(a)+h_j(b)}K_{S_{ij}}(a,b)}{Z_{ij}},
\tag{5}
```

where Z_ij is the normalizer. Its value is the minimum of

```math
\sum_{i<j}\log Z_{ij}-d\sum_{i,a}\nu_i(a)h_i(a).
```

Normalize `r_i(a)=nu_i(a)e^(h_i(a))/sum_b nu_i(b)e^(h_i(b))`. Cancelling
the row normalizers gives exactly (1). This also proves that the infimum
in (1) is attained in the relative interior of the product of local
simplices; infinite divergences exclude their missing-support boundaries.

The needed uniform positivity does not rely on an unproved compactness
claim about m different multipliers. Let `q_i=e^(h_i)`. The aggregate
constraint after cancelling nu_i(a)>0 is

```math
q_i(a)\sum_{j\ne i}
\frac{\sum_b\nu_j(b)q_j(b)K_{S_{ij}}(a,b)}{Z_{ij}}=d.
```

For any two symbols a,a' at i, the two sums inside this expression differ
by a factor between kappa and 1/kappa. Hence `max q_i/min q_i<=1/kappa`.
After normalizing max q_i=1, one has `r_i(a)>=alpha kappa`. Equation (5)
therefore gives the uniform edge-atom lower bound

```math
\pi_{ij}(a,b)\ge\alpha^2\kappa^3.
\tag{6}
```

## 3. Exact-type compiler proves (2)

Let P_* independently sample every edge from its optimizing law (5).
Its aggregate expected row counts equal the desired counts. On T_types,
the likelihood ratio `dP_*/dP_0` is the edge-kernel product times the
constant `exp[-inf Psi_S]`. Thus

```math
Z_S=\exp[\inf\Psi_S]\frac{P_*(T_{\rm types})}{P_0(T_{\rm types})}.
\tag{7}
```

Incidences within a row are independent under P_* because they belong to
different independent edges; their one-site laws need not be identical.
For every symbol the count variance is at most d, so the expected total
number of changes needed to repair all rows is at most `(L/2)m sqrt(d)`.
With probability at least one half, at most `R=ceil(Lm sqrt(d))` changes
suffice. Deterministic rowwise repair alters at most R edges. The uniform
atom bound (6) gives a probability-ratio loss at most
`exp[O_{alpha,K}(R)]`, and each repaired array has at most

```math
\sum_{r\le R}\binom{md}{r}(L-1)^r
=\exp[O_L(m^{3/2}\log m)]
```

preimages. Hence `P_*(T_types)>=exp[-O(m^(3/2)log m)]`, uniformly in S
and all allowed types. Meanwhile independent multinomial type probabilities
give `P_0(T_types)>=exp[-O_L(m log m)]`. Equation (7) proves (2).

## 4. A smooth separated-kernel lemma

**Lemma.** If f is smooth on a neighborhood of `Delta_L x Delta_L`, where
Delta_L is the probability simplex, then there is a finite constant C_f
such that, for all m, real matrices B, and points r_i,s_j in Delta_L,

```math
\left|\sum_{i,j}B_{ij}f(r_i,s_j)\right|\le C_f\beta(B).
\tag{8}
```

Proof: choose a smooth compactly supported extension agreeing with f near
the compact product simplex. After placing its support in the interior of
a cube and periodically extending, its Fourier coefficients c_{k,l} are
absolutely summable. Indeed repeated application of `1-Delta` yields
`|c_{k,l}|<=C_s(1+|k|^2+|l|^2)^(-s)` for s>L. Thus

```math
f(r,s)=\sum_{k,l}c_{k,l}e^{ik\cdot r}e^{il\cdot s},
\qquad\sum_{k,l}|c_{k,l}|<\infty.
```

The bilinear form of B on two complex vectors with coordinate magnitudes
at most one is at most `4 beta(B)`, by decomposing both into their real and
imaginary parts and using the cube extension of beta. Absolute summation
proves (8). The extension and coefficient sum depend only on f, not on m.

This argument may be replaced by a finite-grid approximation if only an
explicit power saving is desired. The smooth lemma avoids a dimension-
dependent loss in the power of m; its constant is not numerically optimized.

## 5. Proof of seed comparison and a seed-free dual

Because K_+ and K_- have entries bounded below, the function

```math
f_{\rm odd}(r,s)=\tfrac12\left[
\log(r^TK_+s)-\log(r^TK_-s)\right]
```

is smooth on a neighborhood of the compact product simplex. It is
symmetric in r,s. With the corresponding half-sum f_even, (1) becomes

```math
\Psi_S(r)=\sum_{i<j}f_{\rm even}(r_i,r_j)
+\sum_{i<j}S_{ij}f_{\rm odd}(r_i,r_j)
+d\sum_iD(\nu_i\Vert r_i).
```

Apply (8) to the hollow symmetric matrix S-T and divide the full ordered
sum by two. The bound is uniform in all r_i, so taking infima proves (3).

In particular define the seed-free functional

```math
\overline\Psi(r)=\frac12\sum_{i<j}
\left[\log(r_i^TK_+r_j)+\log(r_i^TK_-r_j)\right]
+d\sum_iD(\nu_i\Vert r_i).
```

Then

```math
\log Z_S=\inf_r\overline\Psi(r)
+O_{L,\alpha,K}(\beta(S)+m^{3/2}\log m).
\tag{9}
```

The half-sum of logarithms is important: this is NOT the logarithm of the
folded kernel `(K_++K_-)/2`. Removing independent output phases may change
the universal value for asymmetric types. What it cannot do in this fixed-
finite setting is retain the distinction between different bounded-cap
seed constants at leading m² scale.

## 6. Auditing targets and scaling escapes

The delicate points for independent audit are finite-dimensional duality,
the positivity estimate (6), the repair preimage direction in (7), and the
smooth separated-kernel lemma. No conjectural phase statement is used.

The theorem is not uniform as L grows, alpha decreases, or the kernel
condition number grows. Those quantities explicitly control repair and
separated-kernel constants. Thus a possible positive seed-preserving
mechanism must retain growing local complexity, a sufficiently singular
kernel, or genuine cross-row constraints/dependence not represented by
independent permutation types. Such an escape requires its own entropy
and realization accounting; this note asserts no convergence theorem.
