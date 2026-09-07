# Exact affine-equivalence test for the dimension-16 Kerdock probes

Status: proved geometric criterion plus complete finite integer screen.
All 126 five-basis subsets of the proposed nine real dimension-16 Kerdock
bases are affine-net equivalent after a common right rotation and line-row
grouping. This is **not** a classification at growing dimension.

## 1. Criterion invariant under the allowed basis changes

Let U_1,...,U_k be real mutually unbiased orthogonal bases of dimension d=q²,
so every entry of U_i U_j^T is +-1/q. Then the following are equivalent:

1. In some common orthonormal coordinate basis, each U_i consists of q line
   Hadamards supported on q disjoint sets of size q, and any two line supports
   from different families meet in exactly one point.
2. The rows of every U_i can be partitioned into q groups of q such that every
   cross-basis q-by-q tile between two groups has rank one.

The forward implication follows from the one-point intersections. For the
reverse implication, let P_(i,L) be the projection onto the span of the rows
in group L. A cross-group overlap matrix has q² entries of magnitude 1/q,
Frobenius norm one, and rank one. Its singular values are therefore 1,0,...,0.
The two row subspaces intersect in one dimension and their remaining parts
are orthogonal, so their projections commute. Projections inside a family
also commute and sum to I. Thus all line projections commute.

Choose two families. Their q² pairwise intersections are orthogonal rank-one
spaces summing to the ambient space. Their unit vectors give a common point
basis. Every remaining commuting rank-q line projection selects q of these
points, and every cross-family pair intersects in exactly one. Finally, the
left singular vector of any rank-one flat tile has all coordinates +-1/sqrt(q),
so every row's nonzero coefficients in the reconstructed point basis have
that magnitude. The restricted row matrices are exactly line Hadamards.

The criterion only involves U_i U_j^T, so a common right orthogonal rotation
cannot change it. Independent row switches/permutations merely relabel or
resign the partition. Thus this is an exact test against an affine family
in disguise, not a test against one prescribed sparse presentation.

## 2. Complete finite search at q=4

Choose the first basis and enumerate every candidate four-row group containing
row zero: there are binom(15,3)=455. For each other basis, restrict the cross
Hadamard to these four rows and classify its columns by their projective sign
patterns (a pattern and its negative are equivalent).

If a rank-one tiling exists, these classes are exactly its four groups of four.
There cannot be fewer classes because the four restricted rows are independent;
there cannot be more because each of the four desired column groups has rank
one. Thus all other row partitions are forced. Applying the same classification
back through any one other basis also forces the entire first partition.
Checking all pairwise tiles then decides feasibility. No row partition is
missed by this reduction.

## 3. Integer result and explicit common rotations

`tmp/flatify_adversary_2026_09_07_affine_equivalence_screen.py` uses the two
constructive basis generators only for their candidate matrices, and rechecks
every orthogonality and mutual-unbiasedness identity with integers. It recovers
the known five-basis affine construction as a positive control. For all 126
five-element subsets of the nine Kerdock bases it finds valid partitions.
Every tile is tested by the exact sign identity
`T=T[:,0] T[0,:] T[0,0]`, not by a numerical matrix-rank tolerance.

The program additionally reconstructs the actual common rotation for each
subset. Writing U_i=B_i/4, it records an integer matrix W with W^T W=64I,
so the common right rotation is W/8. It verifies every B_i W entry is
0 or +-16, hence every transformed basis entry is 0 or +-1/2; each row has
four nonzeros and all rows in each recovered line group share a support.
The JSON stores all partitions and all these exact common rotations.

As a negative control, the six-basis subset {0,...,5} has no such partition.
This also follows from the general net bound k<=q+1: the centered line-
indicator subspaces have dimension q-1, are mutually orthogonal, and lie in
dimension q²-1. Six families at q=4 are therefore impossible.

## 4. Precise consequence

The dimension-16 five-basis Kerdock diagnostics do not by themselves furnish
a genuinely non-affine escape from the affine-MUB spike mechanism. Their
cross-block matrices are unchanged by the displayed common right rotations,
so they have precisely an affine-net representation.

This finite conclusion does not imply that arbitrary growing-dimensional
five-basis Kerdock families are affine equivalent. It also does not transfer
an asymptotic lower bound solely from one finite example. A separate uniform
geometric argument or a genuine non-affine sequence must be checked on its
own merits. The constructive agent is pursuing such a direct quadratic-phase
argument independently of this finite classification.
