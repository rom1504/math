# Genuine marked old history: exact coherent Walsh root-map bounds

Date: 2026-09-06. Status: new structural lemma, submitted for independent
audit. This note retains the actual marked input and its return. It does
not yet claim the nonlinear feedback energy projection.

## 1. The first marked field and its literal return

Let `B=A/sqrt(n−1)` be a symmetric hollow signing with `||B||op<=L`,
put `Q=B²`, and let `S` be independent signs. Write `m=n−1` and

```math
G=BS,\qquad D_i=S_i h_2(G_i)=S_i(G_i^2-1)/\sqrt2,
\qquad Y=BD,\qquad V=QD.
```

Because `B_ii=0` and `sum_a B_ia²=1`, the input is EXACTLY of Boolean
Walsh degree three:

```math
D_i=\sqrt2\sum_{a<b,\ a,b\ne i}B_{ia}B_{ib}S_iS_aS_b.
```

Its exact covariance is

```math
\operatorname{Cov}(D)=(1-3/m)I+(2/m)Q.                         (1)
```

The diagonal is `1−1/m`; for `i!=j`, the only common marked triples
contain both roots and one further vertex, giving `2 B_ij² Q_ij`.
Thus `Cov(Y)=Q+O_op(1/n)` and `Cov(V)=Q²+O_op(1/n)`.
In particular all individual coordinates of `G,Y,V,QS` have uniformly
bounded second moments. Their Boolean degrees are at most three.

The literal return cannot be treated as a linear-seed coherent field.
It remains pure Walsh degree three and satisfies

```math
\frac{\mathbb E V^{\mathsf T}BV}{2n}
 =\frac{\operatorname{Tr}(B^5)}{2n}+O_L(1/n).                 (2)
```

Hence a first-Boolean-chaos-only energy theorem would discard a
potentially nonzero term. Formula (2) is a scope obstruction, not a
counterexample to a projection that retains `e(V)` literally.

## 2. Exact Boolean gradients

For a vector-valued function `W(S)`, define its gradient matrix by

```math
(J_W)_{ia}=\Delta_a W_i
 =\frac{W_i(S)-W_i(S^{(a)})}{2S_a},
```

where `S^(a)` flips only seed `a`. The entries are independent of `S_a`.
Direct calculation, with no approximation, gives

```math
J_D=\operatorname{diag}(h_2(G))
 +\sqrt2\operatorname{diag}(SG)B
 -\sqrt2\operatorname{diag}(S)(B^{\circ2})\operatorname{diag}(S).
                                                                    (3)
```

Indeed the own derivative is `h2(G_i)`, while for `a!=i` it is
`sqrt2 S_i B_ia(G_i−B_iaS_a)`. Since the absolute row and column sums
of `B^{circ2}` equal one,

```math
\|J_D\|_{op}\le C_L(1+\max_i|G_i|^2),\qquad
J_Y=BJ_D,\quad J_V=QJ_D,\quad J_{QS}=Q.                      (4)
```

The subgaussian maximal bound for `G` yields, for every fixed `p`,
`|| ||J_D||op ||_Lp=O_{p,L}(log(n+1))`; the same holds for `Y,V`.

## 3. Coordinate polynomials preserve polylogarithmic gradient bounds

Let `W^1,...,W^k` be any fixed list selected from `G,Y,V,QS,S`, and
let `P` be a fixed polynomial of total degree `d>=1`. Set
`C_i=P(W_i^1,...,W_i^k)`. The exact finite Taylor expansion gives

```math
\Delta_a C_i
 =\sum_{1\le|\alpha|\le d}
 \frac{(-1)^{|\alpha|+1}2^{|\alpha|-1}}{\alpha!}
 (\partial^\alpha P)(W_i)
 S_a^{|\alpha|-1}\prod_{u=1}^k(\Delta_aW_i^u)^{\alpha_u}.      (5)
```

It follows by writing `W(S^(a))=W(S)−2S_a Delta_a W`. Thus `J_C`
is a finite sum of bounded-row-polynomial diagonals, Hadamard products
of base gradient matrices, and diagonal sign matrices. For arbitrary
square matrices,

```math
\|M_1\circ\cdots\circ M_r\|_{op}
 \le\prod_{u=1}^r\|M_u\|_{op},                              (6)
```

because a Hadamard product is the compression of the tensor product
to the coordinate-diagonal subspace. Therefore (5) bounds the gradient
operator norm by a fixed polynomial of the maximum coordinate sizes
and the base gradient operator norms.

Every base coordinate has degree at most three and uniformly bounded
second moment. Fixed-degree hypercontractivity and the maximal bound
give, for every fixed `p`,

```math
\left\|\max_{i,u}|W_i^u|\right\|_{L^p}
 =O_{p,L,k}((\log(n+1))^{3/2}).
```

Together with (4)--(6), this gives the explicit sufficient estimate

```math
\mathbb E\|J_C\|_{op}^2
 \le C_{P,L,k}(\log(n+1))^{3d-1}.                           (7)
```

For example, a term of derivative order `ell` has maximal factor
`log^{3(d−ell)/2}` and gradient factor `log^ell`; its total exponent
is at most `(3d−1)/2` before squaring. The exact exponent is not
important; every fixed polylogarithm is smaller than any power of `n`.

## 4. Positive Walsh root maps

Discrete Poincare, applied to every deterministic linear combination
of the coordinates, gives the matrix inequality

```math
\operatorname{Cov}(C)\preceq\mathbb E[J_CJ_C^{\mathsf T}].    (8)
```

Write the exact Boolean decomposition `C=sum_q C_q`, including its
deterministic degree-zero row means. Orthogonality of different degrees
gives `Cov(C)=sum_{q>=1} Cov(C_q)`, a sum of positive semidefinite
matrices. Consequently EVERY positive-degree global Walsh root map
`K_q`, in the normalization `K_qK_q^T=Cov(C_q)`, satisfies

```math
\|K_q\|_{op}
 \le C_{P,L,k}(\log(n+1))^{(3d-1)/2}.                       (9)
```

All fixed-root Hilbert norms remain bounded independently of `n`:
the base row moments are bounded by hypercontractivity, and `P` is
fixed. No coherent input-slot collision has been deleted. In
particular, the potentially order-one own-spin contractions inside
`P(QD)` are kept exactly.

The same result holds for cross-covariances, by factorization through
the two root maps. Thus all normalized-energy estimates that require
only bounded root maps and a gain `n^{-c}` remain valid with these
polylogarithmic maps. This does not by itself show that the needed
small-gain contractions exist.

## 5. Next precise gap

For a bounded odd `f(G,Y)`, its first local-Gaussian projection is
`b_0G+b_1Y`. The actual next field decomposes as

```math
Bf(G,Y)=b_0QS+b_1QD+Z,
\qquad Z=B[f(G,Y)-b_0G-b_1Y].
```

The theorem above supplies positive Walsh root-map control for all
fixed polynomial functions of the literal coherent variables in this
decomposition. It does not yet prove local Gaussian separation of `Z`
from `QD`, nor the cross-root mixed-star classification needed for an
energy projection. The planned theorem must retain the full coherent
energy and its bare cross with `Z`; (2) forbids reducing coherent
energy to first Boolean degree.

The numerical falsifier
`computations/continued_feedback_marked_old_tree_projection_2026_09_06.py`
tests the special case `f(Y)=sin(Y)`, using the actual marked field
and actual return throughout. Its numerical agreement is not used in
any proof above.

## 6. A precise fixed-depth corollary

The same argument yields a structural corollary, not a feedback-law
theorem. Start with a fixed number of independent sign-seed vectors.
Allow a fixed finite computation graph whose operations are uniformly
bounded-op deterministic matrix transports and fixed coordinatewise
polynomial maps. Require every input to a matrix transport to have
coordinatewise mean zero. It suffices that each such input be globally
odd under simultaneous seed reversal. Fixed even polynomial auxiliary
fields may be used before a subsequent odd coordinatewise map; they
are not transported without centering.

Then each intermediate field has fixed Boolean degree, coordinate
moments bounded by a fixed power of `log(n+1)`, and gradient operator
moments bounded by a fixed power of `log(n+1)`. Consequently all its
positive-degree exact Walsh root maps have polylogarithmic operator
norms. Constants and exponents may depend on the fixed computation
graph, not on matrix order.

The proof is induction. A polynomial coordinate map uses exactly (5)
and (6), and fixed-degree hypercontractivity bounds its coordinate
moments. A matrix transport has gradient `M J_W`. Its mean is zero by
the transport-input hypothesis, and its coordinate variance is at most
`||M||op² ||Cov(W)||op`, which is polylogarithmic by (8). Fixed-degree
hypercontractivity then bounds its maximal coordinate moments by
another polylogarithm. This supplies the induction assumptions needed
at the next coordinate polynomial. The number and degrees of all
operations are fixed before `n→infinity`.

The centering requirement is essential to this proof. A bounded-op
flat signing can have an apex row whose sum is of order `sqrt(n)`;
transporting the constant vector one then creates a coordinate of
that size, although its input gradient is zero. Thus Poincare alone
does not control transported means.

This corollary removes one root-map estimate from potential fixed
marked-history proofs. It establishes no small proper cuts, source
collision replacement, local Gaussian law, or mixed-star energy
identity for a new circuit. Those remain separate obligations.
