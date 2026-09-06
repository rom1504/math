# The first genuine marked history: local nonlinear-noise separation

Date: 2026-09-06. Status: proof independently reconstructed by the director
and audit agent. This note establishes the local comparison module. The
separate full energy proof is in
`continued_feedback_first_marked_history_energy_projection_2026_09_06.md`;
the local theorem alone does not imply that energy conclusion.

## 1. Frozen setup

Let `B=A/sqrt(n−1)` be a symmetric hollow signing with `||B||op<=L`,
write `Q=B²`, and use independent signs `S`. Set

```math
G=BS,\quad D_i=S_i h_2(G_i),\quad Y=BD,\quad V=QD,
\qquad W_i=(S_i,G_i,Y_i,(QS)_i,V_i).
```

Let `r(g,y)` be a fixed odd polynomial in two standard Gaussian
coordinates, with no first local Hermite chaos. Write

```math
r(g,y)=\sum_{p+q\ge3\ {\rm odd}} r_{pq}h_p(g)h_q(y),
\qquad Z=Br(G,Y),
\qquad T=B\left[\sum_{p,q}r_{pq}^2Q^{\circ(p+q)}\right]B.
```

The proposed local comparison is, for every bounded continuous `Phi`,

```math
\frac1n\sum_i\left|
 \mathbb E\Phi(W_i,Z_i)
 -\mathbb E_{S,N}\Phi(W_i,\sqrt{T_{ii}}N)
 \right|\longrightarrow0,                                  (1)
```

where `N` is independent of the ACTUAL coherent vector `W_i`. The same
holds for polynomial tests, and for bounded functions with allowed
Gaussian-a.e.-continuous factors in `(G_i,Y_i)`, using the established
old Gaussian marginal. No law of `QD` is replaced by a Gaussian.

The covariance `T` and its averaged agreement with the true residual
covariance are supplied by the previously audited old-frame nonlinear
covariance theorem. The new work here is separation from the literal
degree-three return `QD`.

## 2. Exact distinctness preserves a proper-cut bound

Let a fixed-order tensor have any specified left/right flattening.
Imposing equality of one left slot and one right slot is a rectangular
block-diagonal pinching of that flattened matrix. Its operator norm is
at most the original norm: decompose row and column spaces by that
shared label and keep only matched blocks. Imposing inequality is the
identity minus this pinching and costs a factor at most two. A
same-side equality or inequality is a row or column coordinate
projection and costs no factor.

Consequently projection onto pairwise distinct marked labels preserves
every fixed proper-cut operator bound up to a constant depending only
on tensor order. This projection also preserves row Hilbert norms and
the global root-map bound, and commutes with deterministic transport on
the root index. No small Hilbert error from source-diagonal restoration
is assumed here.

## 3. The residual channels have uniformly small proper cuts

The old field `G` has degree-one flat row kernels. The marked field
`Y` has exact degree-three source-Walsh kernels, bounded global root
map, and every proper cut `O_L(n^{-1/2})`. The exact input covariance
is

```math
\operatorname{Cov}(D)=(1-3/(n-1))I+2Q/(n-1).
```

For a local Hermite monomial `h_p(G)h_q(Y)`, the old local Wick
replacement gives the squarefree product of `p` linear `G` branches
and `q` degree-three `Y` branches, in original degree `p+3q`.
Its local row `L²` error is `O(n^{-1/2})` at fixed degrees. The outer
bounded-op transport preserves the averaged `L²` error.

Before source distinctness, the transported product tensor is

```math
K_i^{p,q}=\sum_j B_{ij}\,b_j^{\otimes p}
                                   \otimes y_j^{\otimes q}.
```

The tensor-power cut proof in
`continued_audit_transported_tensor_powers_2026_09_06.md` also applies
to this fixed mixed family. A cut with no straddling old branch is a
product of bounded group root maps with middle `diag(b_i)`. Two or
more straddling degree-three branches give
`sum_j |B_ij| O(n^{-1})=O(n^{-1/2})`. Exactly one straddler is a
block-diagonal partial tensor between bounded whole-branch root maps,
with the remaining root contraction bounded by `||b_i||_2=1`.
There is at least one whole branch because `p+q>=3`.

Thus every proper fixed-root cut is `O(n^{-1/2})`; Section 2 keeps
this bound for the EXACT squarefree source kernel. Global root maps
and row Hilbert norms are bounded. In particular maximum marked-slot
influences are `O(1/n)`.

## 4. A hypergraph contraction lemma retaining coherent collisions

Boolean moment expansion partitions all formal marked positions into
equal-seed blocks of even size. Each primitive exact source-Walsh
tensor has distinct labels internally, so a block meets each primitive
tensor at most once. Blocks meeting three or more tensors must NOT
be discarded merely because they represent non-Gaussian coherent
moments.

The different equality blocks initially carry distinct seed labels.
Remove that inter-block restriction by finite inclusion-exclusion.
Every resulting term is a coarser equality hypergraph; a block that
now meets the same primitive source-Walsh tensor twice gives zero.
All other blocks still have even size at least two. Thus the remaining
argument may sum its hyperedge labels independently without changing
the moment calculation, at the cost of finitely many coefficients
depending only on the fixed total degree.

There is an elementary contraction bound that retains them. Merge two
vertices at once along ALL their common labels. Write `x` for common
labels that also touch a future vertex and must therefore be retained,
and `y` for common labels completed by this merge and summed out. The
merged tensor is

```math
M(x,u,v)=\sum_y A(x,y,u)B(x,y,v).
```

Then

```math
\|M\|_F^2
 =\sum_x\|A_x B_x^{\mathsf T}\|_F^2
 \le\sum_x\|A_x\|_F^2\|B_x\|_F^2
 \le\|A\|_F^2\|B\|_F^2.                                   (2)
```

This covers summing completed labels, retaining shared labels, and
their combination. Repeated vertex merges bound any residual closed
hypergraph by the product of its primitive Hilbert norms; disconnected
components factor. All their original labels occur at least twice, so
no unsummed singleton trace is introduced.

Now suppose `A=K` has all proper cuts at most `epsilon`, and `K`
has MORE formal slots than the other primitive tensor `B`. At least
one slot `u` remains. If `y` is nonempty, each sliced matrix `K_x`
is a row restriction of the proper cut `y | (x,u)`, so
`||K_x||op<=epsilon`. Equation (2) improves to

```math
\|M\|_F\le\varepsilon\|B\|_F.                              (3)
```

If `y` is empty, at least one retained shared label exists; fixing it
uses the maximum influence bound of `K`, again proving (3). This is
why entanglement of the remaining coherent tensors cannot undo the
gain: merge the small-cut noise tensor with one primitive coherent
tensor FIRST, and only then perform the remaining contractive merges.

If the two orders are equal, the same conclusion holds unless they
form an isolated complete pair. When every noise slot is shared but
some labels are retained, slice on a retained label and use its
influence bound. The only exception is the full scalar inner product
of two equal-order tensors with no connection to another vertex.

## 5. The only low-degree complete-pair exception is small

All coherent primitive fields in `W` have original degree at most
three. A residual monomial has degree `p+3q`, with `p+q` odd and at
least three. The sole degree-three possibility is `h_3(G)`; every
other residual degree is at least five.

Let `R_{3,j}` be the exact top Boolean component of `h_3(G_j)`.
Its exact cross with the marked input is

```math
\mathbb E[R_{3,j}D_k]
 =\sqrt3 B_{jk}\left[Q_{jk}^2-
                         \frac{n-2}{(n-1)^2}\right].        (4)
```

For `j=k` both sides vanish. For distinct roots, match the common
triples containing `k`: the coefficient sum is
`2sqrt3 B_jk sum_{a<b} B_jaB_kaB_jbB_kb`, which is (4) because
`sum_a B_ja²B_ka²=(n−2)/(n−1)²`.

The matrix in (4) has absolute row AND column sums `O_L(n^{-1/2})`,
using `sum_k Q_jk²<=L²` and flat entries of `B`. Hence the complete
crosses of `B R_3` with `D,Y=BD,V=QD` all have operator norm
`O_L(n^{-1/2})`. A complete pair with a degree-one coherent field is
impossible. This handles every exceptional isolated noise/coherent
pair from Section 4.

## 6. Polynomial local independence

Expand a fixed polynomial moment in the coherent primitive fields and
the finitely many exact transported residual channels. The primitive
noise tensors have uniformly small proper cuts. Any connected moment
component meeting noise and coherent vertices gains `O(n^{-1/2})`
by the first-merge argument and Section 5. Any component containing
more than two noise vertices has a proper or retained-label overlap
and the same gain. The only surviving noise components are isolated
complete pairs. Their contractions are their actual covariance values.

Therefore the joint polynomial moments factor into the ACTUAL Boolean
coherent moments and independent Gaussian noise pairings, with a
uniformly vanishing fixed-degree error. Different original input
degrees remain exactly orthogonal. If several residual monomials
have the same original degree, retain their actual within-degree
covariance rather than declaring them orthogonal without proof.

For the total noise, let `v_i` first denote the exact squarefree
variance. The already established nonlinear covariance comparison
gives `n^{-1}sum_i|v_i−T_ii|→0`. This changes averaged fixed polynomial
Gaussian tests by `o(1)` because the variances are uniformly bounded.
The local-Wick source errors are removed after transport using their
averaged `L²` bound and fixed-degree hypercontractivity. At fixed
polynomial degree their row moments remain bounded: the summed
source-error variances are `O(1)`, and `||b_i||_2=1` bounds each
transported error variance by that sum. Higher moments follow by
hypercontractivity. This proves the polynomial version of (1).

## 7. The literal return has a uniform exponential moment

A general degree-three chaos bound alone would give only sub-Weibull
tails and would not justify moment determinacy. The marked structure
gives a stronger bound. For any deterministic row weights `w` of
bounded Euclidean norm, put `V_w=sum_j w_j D_j`. Randomly color the
original seed indices independently and uniformly into three classes
`C_0,C_1,C_2`. The exact tetrahedral expansion gives

```math
V_w=\frac{27}{\sqrt2}\mathbb E_{\rm color}
 \sum_{j\in C_0}w_jS_j
   \left(\sum_{a\in C_1}B_{ja}S_a\right)
   \left(\sum_{b\in C_2}B_{jb}S_b\right).                  (5)
```

Every ordered triple of distinct labels receives the specified three
colors with probability `1/27`. Condition on the colors and the
`C_0` seeds. The remaining sum is a bilinear form of independent
Rademacher vectors, with matrix

```math
M=B_{C_0,C_1}^{\mathsf T}\operatorname{diag}(w_{C_0}S_{C_0})
 B_{C_0,C_2},\qquad \|M\|_F\le L\|w\|_2.
```

Two Khintchine inequalities, with Minkowski for the intermediate
Euclidean norm, give `||X^T M Y||_Lp<=(p−1)||M||F` for `p>=2`.
Jensen in the coloring average therefore proves

```math
\|V_w\|_{L^p}\le (27/\sqrt2)L\|w\|_2(p-1).                 (6)
```

This implies a uniform exponential moment near zero. Taking `w` to
be rows of `B` or `Q` proves this for `Y_i` and `V_i`. The other
coherent primitives are bounded or subgaussian. Their fixed-dimensional
joint law consequently has a uniform exponential moment in the sum
of absolute coordinate values.

## 8. Bounded local closure and limitation

The exponential-moment bound makes each subsequential coherent limit
moment-determinate. The squarefree noise has asymptotically Gaussian
moments and uniformly bounded fixed moments. Any subsequential joint
limit therefore has an exponential moment near zero as well: its
coherent marginal does, its Gaussian noise marginal does, and
Cauchy--Schwarz controls their joint exponential moment without first
assuming independence. The factored moments in Section 6 force this
joint limit to be the product of its literal coherent marginal and
the comparison Gaussian. Compactness and a worst-root subsequence
argument give the uniform bounded-continuous comparison for squarefree
channels. The averaged raw replacement and variance comparison give
(1) as stated.

For a bounded odd Gaussian-a.e.-continuous `f(g,y)` (or an explicitly
stipulated ordered actual-marginal `L²` realization), subtract its first local Gaussian
projection and approximate the residual by fixed Hermite polynomials.
The old Gaussian marginal supplies the averaged `L²` approximation,
the fixed operator cap transports it, and independent Gaussian channel
coupling controls the comparison variances. Thus the bounded-continuous
version extends in the ordered limit, including fixed Lipschitz
feedback functions and Gaussian-a.e.-continuous bounded old masks.

This proves only a LOCAL separation from the literal return. An actual
self-energy theorem still needs the cross-root mixed-star classification
for the same marked coherent tensors. A local law alone cannot
identify `B[H(G,Y) psi(Bf(G,Y))]` or its quadratic energy.
