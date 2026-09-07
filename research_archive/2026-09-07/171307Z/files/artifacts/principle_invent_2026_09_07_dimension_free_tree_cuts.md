# Dimension-free injective tree cuts under entry delocalization

Date: 2026-09-07. Status: **proved cut and collision lemmas; not a full
extension of the marked lower theorem**.

Let B be a real symmetric hollow matrix, with row square sums at most r,
operator norm at most L, and b=max|B_ij|. Constants below depend only on
the fixed tree, r,L. For a tree T, its coefficient tensor is the product
of B entries on its edges, one index per vertex. A marked/output tree
has an external root of degree one.

## 1. Global cuts need no ambient-dimension entry bound

For any nontrivial bipartition of the vertices into row and column
indices, the crossing edges give a tensor product of B matrices. Its
norm is at most max(1,L)^|E(T)|. Copying the label of a crossing-incident
vertex to its crossing-edge slots is an isometry.

On either side of the cut, consider the internal-edge forest. Every
component contains a crossing-incident vertex, since T is connected.
Regard the labels of these incident vertices as boundary labels. The
map that expands a fixed boundary labeling to all labels on that side,
weighted by the product of internal B entries, has disjoint supports
for different boundary labelings. Its operator norm is thus the largest
square root of a squared column sum.

That squared column sum is the sum of a product of P=B circ B entries
over the nonboundary vertices. Choose one boundary root in each internal
component and successively sum tree leaves. Additional boundary labels
are fixed instead of summed, which can only decrease this nonnegative
sum. The result is bounded by max(1,r)^|E(T)|.

Composing these expansions, copy maps, and the crossing-edge tensor
product proves a global flattening bound C_T,r,L independent of n and b.
This replaces the coarser sqrt(n)-per-isolated-vertex and maxentry-per-
internal-edge argument.

## 2. Fixing an external root gives a factor b

Fix the degree-one external root to label i, and make a proper cut of
the remaining marked vertices. Assign the external root to the same
side as its neighbor. Its incident edge is then internal. The internal
component containing this fixed root has a crossing-boundary vertex,
at some positive graph distance ell.

For fixed label j at one such boundary vertex, summing the squared
weights along the root-to-boundary path gives (P^ell)_ij. Since P is
nonnegative, its row sums are at most r, and maxentry(P)<=b^2,

    (P^ell)_ij <= b^2 r^(ell-1),       ell>=1.

The other branches and components contribute only fixed powers of r.
Thus the expansion map on this side has norm O_T,r(b), while the other
maps have bounded norm. Every proper marked-slot flattening at a fixed
external root is consequently O_T,r,L(b).

## 3. Injectivity masks preserve the cut bounds directly

No large global Hilbert error is needed to pass to injective tensors.
For one forbidden equality of two labels:

* If the two slots are on the same side of a cut, its equality indicator
  is an orthogonal row or column projection.
* If they lie on opposite sides, its equality indicator acts as
  M -> sum_v P_v M Q_v, where (P_v) and (Q_v) are orthogonal coordinate
  partitions. This rectangular block-diagonal compression has operator
  norm at most one.

Removing an equality therefore has multiplier norm at most two. There
are only a fixed number of slot pairs. Their entrywise masks commute,
so their product changes every cut norm by at most 2^(number of pairs).
Equalities with a fixed external label are ordinary row/column masks.
The global O(1) and fixed-root O(b) bounds thus hold for the fully
injective tree tensor itself.

## 4. Weighted collision mass is small at each root

When P has exact row sums one, the squared unrestricted tensor rooted
at i is the law of a tree-indexed Markov chain with transition P. For
two distinct tree vertices at distance ell, the probability that their
labels coincide, conditional on the label at their last common ancestor,
is (P^ell)_vv. Symmetry of P gives the same formula for vertices in
different descendant branches. It is at most b^2.

The union bound over the fixed number of vertex pairs shows that
collision deletion has squared Hilbert norm O_T(b^2) at every fixed
root. Its global squared Hilbert norm is O_T(n b^2), which is small
after division by n if b->0, but need not be O(1).

Uniform approximate row sums 1+o(1) give the same bounds up to bounded
fixed powers and o(1) normalization changes. This proves exactly the
collision estimate needed for normalized nuclear covariance transfers,
without pretending that its global Hilbert norm is dimension-free.

## 5. What remains open in this attempted strengthening

These lemmas extend the generic cut mechanism from max|B_ij|=O(n^-1/2)
to max|B_ij|=o(1), at fixed normalized operator norm and row norms.
They do NOT by themselves establish the finite old Gaussian frame or
its covariance-operator Gram identity for all such matrices.

The dense weighted proof uses bounded W=sqrt(n)B and normalized
single-edge diagram cancellation. That argument is not available
unchanged when sqrt(n)b grows. Nonisomorphic whole-tree overlaps and
old-input Gram off-diagonals must instead be bounded by a genuinely
dimension-free graph/tensor inequality. Generic small proper cuts are
not a substitute for proving these full-overlap covariances.

For the first nontrivial old input, the covariance does have a helpful
form: its off-diagonal is a fixed multiple of

    (B circ B) circ (B^2).

Its absolute row sums are at most b L, by Cauchy--Schwarz and
sum_j B_ij^4<=b^2 sum_j B_ij^2. This verifies this one frame coordinate,
not the arbitrary-tree claim.

No dimension-free extension of the numerical lower bound is asserted
until those full-overlap modules are proved or a counterexample is found.

## 6. Preserved exact full-overlap graph gap at the priority change

An off-diagonal input Gram pattern is the union, with multiplicity, of
two spanning trees on v vertices. Every cut therefore has at least two
edge occurrences. The two output roots i,j have odd total degree, and
every other vertex has even degree. Every i-versus-j cut is odd and
has at least three edges, giving three edge-disjoint root paths.

There is also always a free degree-two vertex. Otherwise all v-2 free
vertices would have degree at least four and the two roots at least
three, totaling at least 4v-2 degrees, whereas the two trees have only
4v-4 degree occurrences. Such a free vertex is a leaf in BOTH original
input trees. A common parent allows an exact doubled-leaf summation;
distinct parents replace the two leaf edges by Q_ab.

The latter operation produces a Q-decorated pair of pruned trees and
destroys the simple parity-degree induction. The existing prescribed-
boundary Mingo--Speicher graph lemma supplies O(1), not automatically
the needed o(1), for the resulting root operator. Three root paths
alone are insufficient: three length-two paths give Q^(circ 3), which
need not be small. The two-spanning-tree edge count adds a real condition
but has not yet been converted into the required quantitative bound.

For the first nontrivial Q-decorated theta model, put
z_ij(u)=B_iu B_ju and P=B circ B. Its kernel is

    K_ij=Q_ij z_ij^T P z_ij.

The stochastic contraction of P gives |z^T P z|<=||z||^2;
max_j ||z_ij||^2<=b^2 and sum_j ||z_ij||^2=1. Thus row Cauchy--Schwarz
gives sum_j |K_ij|<=L b. This model passes. It is not a proof for all
decorated two-tree graphs.

Main research priority changed at the director's request to a concrete
lower/upper Gaussian-control duality test. The full delocalized lower
bound remains OPEN, while the bounded-amplitude result is separately
audited and proved.
