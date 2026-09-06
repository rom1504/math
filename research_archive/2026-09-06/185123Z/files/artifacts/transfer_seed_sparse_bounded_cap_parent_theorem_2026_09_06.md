# Sparse random restrictions of every bounded-cap parent sequence

Date: 2026-09-06. Seed-transfer track. This is a statement about ACTUAL
signings with an actual cap bound; no minimizing-seed or Hadamard
assumption is made in the final theorem.

## Main theorem

Let `A_D` be any sequence of symmetric hollow sign matrices with
`Q(A_D)=O(D^(3/2))`. Let `n=n(D)` tend to infinity with `n/D->0`, and
let `T` be a uniform subset of exact size `n`. Then, for every `eta>0`,

```math
\Pr\left\{\frac{Q((A_D)_T)}{n^{3/2}}
                  \ge\frac2\pi-\eta\right\}\longrightarrow1. (1)
```

Thus uniform random principal restriction at vanishing retention does
not preserve the strict-subhalf caps of the known good parent family,
or of any asymptotically minimizing family. The obstruction concerns
TYPICAL selectors. It does not rule out exceptional deterministic
selectors, give a lower bound `2/pi` on `M_n`, or settle convergence.

The proof has three distinct ingredients: a local-moment theorem for
bounded-operator parents; the finite Gaussian covariance witness; and
the already-proved principal spectral deletion. Their scopes are kept
separate below.

## 1. The bounded-operator local-moment theorem

Assume first that `||A_D||op<=C sqrt(D)` for a fixed `C`. For uniform
`T` of size `n<=D`, set `L=(A_D)_T/sqrt(n)` and `p=n/D`.
For every fixed `k`, with semicircle moments
`m_(2j)=Catalan(j), m_(2j+1)=0`,

```math
\mathbb E\frac1n\sum_i\big((L^k)_{ii}-m_k\big)^2
=O_{k,C}(p+1/n+D^{-1/2}).                                 (2)
```

In particular (2) is uniform in simultaneous `n->infinity`, `p->0`.
No spectral assumption on the unregularized parents is used yet.

### 1.1 Exact imported graph bound and its hypotheses

The primary result used here is [Mingo--Speicher, *Sharp Bounds for
Sums Associated to Graphs of Matrices*, Theorem 6, pp. 6--7](https://arxiv.org/pdf/0909.4277).
It bounds an unrestricted matrix-entry graph sum by a power of the
matrix size times the product of the edge operator norms. The exponent
is determined by the leaves of the forest of two-edge-connected
components. A component with no cut edge contributes exactly one to
that exponent. The theorem explicitly permits loops and multiple
edges; changing an edge orientation transposes its matrix and does not
alter the operator norm.

Every connected Eulerian graph has no cut edge: summing degrees on one
side of a hypothetical cut edge would have odd parity. Therefore, for
an Eulerian graph `F` with `v` vertices, `e` edges, and `c` components
(isolated vertices included), the theorem specializes to

```math
\left|\sum_{z:V(F)\to[D]}\prod_{uv\in E(F)}B_{z_u,z_v}\right|
\le D^c\|B\|_{op}^{e}.                                    (3)
```

The same primary input-output graph machinery was already used in
`continued_director_boundary_graph_feedback_cuts_2026_09_06.md`; the
new use here is its unrestricted Eulerian graph-sum corollary.

### 1.2 Full-sign completion, parity, and the zero-excess case

Set `B=A_D+I`. This is a full symmetric SIGN matrix, not the normalized
operator. Its norm is at most `(C+1)sqrt(D)`, and hollowing a principal
restriction of `B` returns exactly the target restriction of `A_D`.
Using this completion is essential: cancellation of even powers is
exact for full signs, but would fail on coincident colors at a hollow
diagonal.

For a graph `F`, put

```math
t_F(B)=D^{-v}\sum_z\prod_{uv\in E(F)}B_{z_u,z_v},\qquad
b(F)=v-c-e/2.
```

Equation (3) gives, for every Eulerian `F`,

```math
|t_F(B)|\le(C+1)^eD^{-b(F)}.                               (4)
```

This replaces the exact Hadamard degree-two elimination. It is NOT
obtained by replacing a Kronecker identity with an inequality inside
that previous elimination proof.

A second estimate handles the missing case `b=0`. If `F` has any
nonloop edge `uv`, condition on all colors except those at `u,v`.
The remaining product is a fixed sign times `B_ab f(a)g(b)`, where
`f,g` are sign functions. Hence

```math
|t_F(B)|\le\beta(B)/D^2
\le\|B\|_{op}/D\le(C+1)D^{-1/2}.                           (5)
```

For the last inequality `beta(B)<=D||B||op` suffices; no original cap
bound or Grothendieck inequality is needed in this bounded-op step.

### 1.3 Diagram and injectivity proof

Use the same closed-walk notation as in
`transfer_seed_sparse_hadamard_local_moments_2026_09_06.md`: a connected
walk of length `k` has `v` vertices and exponent `a=v-1-k/2`; reduce
edge multiplicities modulo two to its Eulerian parity graph `F`.
The component-connection counting argument gives `a<=b(F)`.

Thus when `a>0`, (4) bounds the normalized contribution by
`(C+1)^k p^a`. When `a<0`, its absolute value is at most `n^a`.
If `a=0` and `F` is nonempty, then `F` must contain a nonloop edge:
a nonempty union of loops and isolated vertices has negative `b`.
Equation (5) therefore makes this remaining contribution
`O_C(D^(-1/2))`.

The only leading pattern is consequently an empty parity graph with
`a=0`: a tree with every edge used exactly twice. It has sign product
one. Its contour walks give the Catalan moments.

Exact-size ambient injectivity is handled by the same partition-Mobius
identity as before, not by approximating the entire sample by distinct
independent words. A partition identifying `h>=1` colors contributes

```math
n^aD^{-h}t_{F/\pi}(B)
=p^h n^{a-h}t_{F/\pi}(B).                                 (6)
```

The quotient is another connected closed-walk multigraph, now allowing
loops. Its exponent is `a-h`, and (4) plus the trivial bound give
`|n^(a-h)t_(F/pi)|<=max(1,(C+1)^k)`. Every such correction is thus
`O_{k,C}(p)`. All sample-label and denominator corrections are
`O_k(1/n+1/D)`.

For a pair of rooted closed walks, the leading doubled tree can use
an edge in only one of the two walks: each closed walk separately
uses every tree edge an even number of times. Their connected edge
sets therefore meet only at the root. The leading count factors into
the product of their Catalan counts. As in the Hadamard proof, the
rooted second moment proves (2). All errors at even combined length
are `O_{k,C}(p+1/n+D^(-1/2))`.

## 2. Fixed-degree Gaussian rounding reaches two over pi

For completeness, use the polynomial hierarchy from
`transfer_adversary_sparse_hadamard_two_over_pi_2026_09_06.md`.
Let `U_j(x/2)` be the variance-one semicircle orthonormal polynomials,
specified by `U_0=1,U_1=x,U_(j+1)=x U_j-U_(j-1)`. For a FIXED `q`, set

```math
P_q(x)=\sum_{j=0}^qU_j(x/2),\qquad R=I+P_q(L)^2.
```

The semicircle moments give `tr(R)/n->q+2` and
`tr(LR)/n->2q`. Equation (2) supplies diagonal concentration, while
the required higher fixed moments bound `tr(R^2)/n`. Since `R_ii>=1`,
the finite covariance inequality from
`transfer_seed_finite_spectral_gaussian_cap_2026_09_06.md` gives

```math
Q((A_D)_T)/n^{3/2}\ge\frac{2q}{\pi(q+2)}-o_P(1).            (7)
```

First let `D` tend to infinity with FIXED `q,C`, and only then let
`q` increase. This proves (1) for bounded-operator parents. There is
no assertion of a uniform growing-degree polynomial approximation.

## 3. Removing the operator assumption using actual cap information

The principal spectral-deletion theorem is proved in Section 1 of
`resumed_bound_audit_minimal_proof_2026_09_06.md`, independently of
the numerical universal lower bound. For every fixed `epsilon>0`,
it supplies a principal core `R_D` of size at least `(1-epsilon)D`
such that

```math
\|(A_D)_{R_D}\|_{op}
\le\frac{4K_G Q(A_D)}{\epsilon D}
=O_\epsilon(\sqrt D),\qquad
K_G=\frac\pi{2\operatorname{arsinh}(1)}.                    (8)
```

Its inputs are polarization `beta(A)<=4Q(A)`, a proved Grothendieck
diagonal majorant, and deletion of the large diagonal entries. It does
not assume original near-minimizers have bounded operator norm.

For a uniform `n`-subset `T`, let `m=|T intersect R_D|`. Hypergeometric
variance is at most `n/4`, so
`m/n>=1-epsilon-o_P(1)`. Conditional on `m`, the intersection is a
uniform `m`-subset of `R_D`. Also `m->infinity` in probability and
`m/|R_D|->0`. The quantitative fixed-degree moment bounds are uniform
over the likely range of `m`, so Section 2 applies to these cores.

Principal monotonicity is exact: averaging all omitted unbiased spins
shows `Q(A_T)>=Q(A_(T intersect R_D))`. Therefore, for fixed
`epsilon` and fixed polynomial degree `q`,

```math
\frac{Q((A_D)_T)}{n^{3/2}}
\ge\left(\frac{2q}{\pi(q+2)}-o_P(1)\right)
                       (1-\epsilon-o_P(1))^{3/2}.          (9)
```

Choose `epsilon` small and `q` large for a desired fixed error, then
take the order limit. This proves the main theorem (1) for every
bounded-cap parent sequence.

## Remaining scope

The theorem says a uniformly selected small fraction of a good large
signing is typically substantially WORSE after its own normalization.
It is compatible with exceptionally good principal restrictions and
with the all-order strict-subhalf upper construction. It supplies no
method for finding those exceptional selectors and no relation forcing
the original liminf and limsup to coincide.
