# Boundary-covered graph kernels and high-degree marked returns

Date: 2026-09-06. Director proof, independently reconstructed by two agents.
This is an exact finite-polynomial estimate. The bounded-response passage,
Gaussian regression, and useful innovation size are separate obligations.

## 1. A matrix-kernel lemma with prescribed boundary sides

Let G be a finite connected undirected multigraph, with loops permitted.
Each edge e carries a real square matrix M_e on a common n-point set.
Transposing a matrix changes its orientation but not its operator norm.
Let a nonempty set of vertices be boundary vertices. For fixed boundary
labels, sum the product of edge entries over all other vertex labels;
the resulting tensor is K_G. Multiple edges are retained, not merged.

Assume every leaf of the forest obtained by contracting the two-edge-
connected components contains a boundary vertex. A trivial one-vertex
forest counts as a leaf. For every partition of the boundary into two
nonempty sets L,R, the corresponding matrix flattening satisfies

```math
 \|K_G\|_{L\mid R,op}\le\prod_e\|M_e\|_{op}.             (1)
```

If an additional factor w(label(v)) is inserted at any vertex v, the
right side is multiplied by ||w||_infinity. Several such factors can
be inserted, with the product of their bounds. There is no factor n
for a closed cycle attached to an internal vertex.

### Proof of the boundary/orientation point

Attach one external identity-edge port at each boundary vertex. Mark
ports in L as inputs and those in R as outputs. Introduce temporary
vertices s,t, edges from s to all input ports, edges from all output
ports to t, and an edge st. This augmentation is used ONLY to choose
an orientation: s,t have no matrix space and do not identify different
boundary labels.

The augmented graph has no bridge. For an original bridge, each side
contains a boundary vertex by the leaf hypothesis; the temporary hubs
and st provide an alternative connection. Each port edge and each
hub edge also lies on a cycle, since both boundary colors are present.
Non-bridge edges remain non-bridges.

The input-output modification argument of Mingo--Speicher, Lemma15,
now orients this graph acyclically from s to t after splitting internal
vertices with identity edges. One detail matters: s and t must not be
split. Both augmented G-s and G-t are connected. When adding an unused
edge incident to s, its return path can be chosen to reach the current
graph through a path to t avoiding s; analogously at t. Thus no cycle
ear forces a split of either hub. The external ports have degree two
and also need no split. Delete both hubs and all temporary hub/st
edges. The remaining network has the prescribed input and output
ports, and every internal vertex has an incoming and outgoing edge.

At each internal vertex use the equality map

```math
 \sum_j |e_j^{\otimes a}\rangle\langle e_j^{\otimes b}|,
 \qquad a,b\ge1.
```

It is a partial isometry. Edge maps, tensor products, and these maps
in a topological order realize K_G exactly. Their norms multiply to
(1). A vertex weight replaces one equality map by its weighted version,
whose norm is ||w||_infinity. If a vertex was split, put the weight at
one copy; identity constraints force all its copies to have equal labels.

Primary source for the input-output norm bound and identity-edge ear
modification: [Mingo and Speicher, Sharp bounds for sums associated to
graphs of matrices, Theorem11 and Lemma15](https://arxiv.org/abs/0909.4277).
The prescribed two-sided boundary formulation above includes the hub-
preservation argument; a generic partial transpose of a previously
bounded graph operator would NOT prove it.

## 2. Exact diagrams for the marked primitives

Let B be symmetric and hollow, with unit squared row sums, operator norm
at most L, and max entry magnitude at most epsilon. Set Q=B^2. On the
uniform Boolean cube let

```math
 G=BS,\quad D_j=S_j h_2(G_j),\quad Y=BD,
 \quad W=(S,G,Y,QS,QD),\qquad h_2(u)=(u^2-1)/\sqrt2.
```

For any fixed polynomial P in the five coordinates of W, take the
EXACT Walsh expansion of P(W_j). Its primitive diagrams are:

- S_j: one seed mark at the root j, no edge;
- G_j or (QS)_j: one B or Q edge from j to a marked seed;
- Y_j or (QD)_j: one B or Q edge from j to k, two B edges from k
  to a,b, and seed marks k,a,b, required to be distinct.

The last assertion is exact, since

```math
 D_k=\sqrt2\sum_{a<b}B_{ka}B_{kb}S_kS_aS_b.
```

Unit row norm cancels the constant term of h2, and hollowness removes
the seed k from a,b. Fixed factorial and Hermite constants can be
absorbed into constants depending on the polynomial degree.

Multiply a fixed finite number of these diagrams. Partition the seed
marks according to equality of their labels. Blocks of odd cardinality
are the surviving free Walsh labels; even blocks are summed internally.
This is the ordinary exact Boolean identity S_a^2=1, not a Gaussian
Wick approximation. Internal centers coincide with their seed mark,
so every vertex other than the common root has graph degree of the
same parity as its number of seed marks. Thus every internal non-root
vertex has even degree, and each free non-root vertex has odd degree.

Distinctness of different free labels is imposed at the end. Distinctness
between internal labels, or between internal and free labels, is handled
by finite inclusion-exclusion. If a coarsening transitively merges two
declared free vertices, its term is killed by the final free-distinct
projection and is discarded. Every other coarsening preserves the parity
description and the number of free vertices. The matrices are still B,Q,
with multiple edges and loops retained.

The root is not silently excluded from seed labels. If S_j is one of
the factors, its mark participates in the seed equality partition.
If the root carries no seed mark, its summation label remains unrestricted
and can equal a free or internal seed label. In particular the diagonal
contribution Q_jj S_j to (QS)_j is retained. Equivalently one can refine
the bookkeeping by root equalities without deleting any such terms.

## 3. The small-cut theorem

Write K_{i,q} for the exact degree-q Walsh coefficient tensor of

```math
 V_i=\sum_j B_{ij}P(W_j).
```

For every fixed q>3 and every nontrivial partition of its q seed slots,

```math
 \max_i\|K_{i,q}\|_{cut,op}
       \le C_{P,L}\,\epsilon.                           (2)
```

The coefficient convention (ordered squarefree tuples or subsets with
symmetrization) changes only C_{P,L}. This statement concerns proper
cuts, not the Euclidean norm of the entire coefficient vector.
The same statement holds for root-dependent polynomials P_j of uniformly
bounded degree and coefficients: each coefficient is another bounded
weight at the common root. This includes the deterministic root-dependent
coefficients used at a fixed response-approximation stage.

### Why every diagram has enough boundary

Consider a connected equality diagram before the outer row B_ij is
inserted as a weight at root j. A bridge side not containing j includes
at most one original primitive center: each such center has its own
edge directly to j, so two such centers would supply two crossing edges.
If it has one center, only that primitive's at most three seed marks
can occur there. If it has no center, it contains only a single original
leaf mark. Identified copies retain multiple root edges and cannot
evade this count. Hence such a side has at most three distinct free
Walsh vertices.

A bridge-forest leaf component away from the root must contain a free
vertex. Otherwise all its vertices would have even degree but its
single exiting edge would make their total degree odd. The only
possible exceptional leaf component is the one containing j. If that
component has no free vertex, all q free vertices lie on the other
side of its unique bridge; the preceding count gives q<=3. For q>3
this is impossible. The trivial-forest case contains the q free
vertices automatically. Therefore every bridge-forest leaf is covered
by a free boundary vertex.

Apply (1) for any proper partition of the free vertices, with the
root weight w_j=B_ij. Its norm is at most epsilon. Every edge matrix
has norm at most max(L,L^2), and the number of edges and equality
patterns is fixed by P. This proves (2) for each internal-distinctness
diagram and then their finite sum.

Finally exclude equal free coordinates. Equality on one side of a
matrix cut is a coordinate projection. Equality across its two sides
is block pinching, a contraction in operator norm. Its complement
has norm at most two. Applying finitely many such exclusion maps only
changes C_{P,L}; it does not sum over free labels or create a dimension
factor. This completes the exact squarefree estimate.

## 4. What it would remove in the current feedback proof

Take P=c0 at a fixed polynomial approximation stage. Equation(2)
controls every original Walsh degree above three in Bc0(W), without
pretending QD has small proper cuts. If a new-noise degree p>=5 is
fully contracted into a degree q>p of Bc0, that contraction is a
proper cut of the RIGHT tensor. Its squared norm is bounded by
C epsilon^2 times the squared norm of the left tensor. Thus bounded
averaged left variance supplies the required averaged o(1) contraction.

This is a finite-stage consequence, not a claim that a bounded response
has a fixed finite polynomial degree. Limits of approximations must
remain ordered. The other retained return BD_aB r, equal-degree
regression, low-variance truncation, and positive innovation size are
not replaced by (2). No new universal constant or original convergence
theorem follows from this estimate alone.

## 5. The exact structural generalization

The same proof applies to a fixed finite collection of rooted matrix
graph primitives with at most d seed marks each, provided:

1. each primitive is connected to its distinguished root;
2. at every non-root vertex its graph degree has the parity of its
   number of seed marks;
3. all edge matrices have uniformly bounded operator norms;
4. multiplying primitives retains their separate edges, even when
   equality of labels identifies their endpoints.

For a bounded-degree polynomial in these primitives, followed by a
matrix with maximum entry epsilon, every proper cut of its exact
Walsh degree q>d is O(epsilon). The final matrix need not be symmetric
for this assertion. Constants depend on the fixed primitive graphs,
polynomial, and edge norm bounds, not on n.

Indeed, if a bridge side away from the common root contained a seed
mark from two different primitive copies, the two original root-to-mark
paths, made of disjoint copies of edges, would cross that cut in two
different edges. Thus it contains marks from at most one primitive,
at most d in total. The same even-degree and boundary argument proves
the result. This is a graph-parity hypothesis, not an assertion about
arbitrary bounded-op polynomial arrays: an arbitrary pre-mixed high-degree
chaos can retain non-Gaussian proper cuts after a second transport.

## 6. Global diagram cuts and the even-degree sharpening

For any individual equality diagram of positive Walsh degree, declare
the source root itself to be an additional boundary port. The root
bridge-leaf is now covered automatically; a non-root uncovered leaf
is again impossible by parity. Applying (1) proves O_{P,L}(1) for
EVERY proper global cut of that diagram, including cuts with the root
on either side. If the root equals a seed vertex, retain separate
identity ports at the same vertex; their equality gates are contractions.
Summing finitely many diagrams gives the same all-global estimate for
the exact polynomial's positive-degree tensors. Degree zero is excluded:
its vector of root means can have Euclidean norm of order sqrt(n).

There is also no exceptional positive EVEN degree in (2). A globally
even source monomial has even total seed-mark count, so the parity
identity at the root agrees with its own mark count as well. Every
internal vertex, including an unmarked root, then has even degree.
Thus every bridge leaf contains a free mark without needing q>3.
All positive even degrees (in particular q=2) have O(epsilon) proper
cuts after flat transport.

Degree three is genuinely exceptional. On a conference sequence,
Q=I and P(W)=Y gives BY=D=S h2(G). The proper cut separating its
own-spin mark from the two remaining seed slots has a norm bounded
away from zero. Hence the cubic threshold is not merely a weakness
of the graph proof. Degree one has no nontrivial proper seed cut.

[Independent proof and fixed-stage next-return consequence](continued_audit_boundary_graph_and_next_return_2026_09_06.md).
