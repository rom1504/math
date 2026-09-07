# Bipartite extension audit of the marked Gaussian lower certificate

Verdict: the inspected old-frame, tensor/forest, exceptional-return,
odd-Schur, and cap-conversion mechanisms survive for COMPLETE BIPARTITE
full-sign support. No bipartite parity obstruction was found. This note
does not assert an extension to arbitrary small-entry row-normalized
weighted matrices.

Let C be an arbitrary m by m full-sign matrix and set

    A = [[0,C],[C^T,0]],  N=2m,  B=A/sqrt(m),  Q=B^2.

Then B is hollow symmetric, B has unit squared row norm, max|B_ij|=m^-1/2,
and Q is a PSD correlation matrix (block diagonal is allowed). Denote
L(C)=max_(x,y Boolean)|x^T C y|. Exactly

    Qcap(A)=L(C),   beta(A)=max_(x,y)|x^T A y|=2L(C).   (1)

The candidate consequence is

    liminf_m min_C L(C)/m^(3/2) >= 2 c_*,
    c_*=.43332211166408075341581292889757934655863464803369... . (2)

The scalar numerical certificate is unchanged; the proof obligations here
are its matrix/tree interface and the square-core reduction.

## 1. What changes in parity-edge diagram arguments

Write A_ij=0 on same-shore pairs and A_ij in {+-1} on opposite shores.
For odd k, A_ij^k=A_ij. For positive even k, A_ij^k is the support
indicator of an opposite-shore pair. Thus one must NOT simply erase
the even powers as if all off-diagonal entries were signs.

Instead keep the support constraints. A connected underlying diagram
which is not bipartite contributes zero. If it is bipartite, fixing the
shore of one root uniquely determines the shore of every other vertex;
additional fixed roots either satisfy this coloring or make the term
zero. Each free vertex then has m possible labels rather than N. At a
fixed number of vertices, injectivity exclusions cost O(1/m).

If an odd-multiplicity edge with two free endpoints survives, fix all
other labels. All remaining factors separate into bounded unary weights
at the two endpoints; the only joint factor is A_uv. The sum is bounded
by beta(A), including the shore indicators and forbidden-label exclusions
in the unary weights. Dividing by m^2 gives the same o(1) discrepancy
gain under L(C)=O(m^(3/2)). Graphs in doubled-Frobenius estimates work
the same way, with all root labels summed.

If the parity graph is empty at leading label count, the original edge
count/connectivity argument still forces a doubled tree and whole-copy
pairing. Its normalized squared-edge weight is exactly one before
injectivity exclusions: remove leaves using sum_j B_ij^2=1. The rooted
automorphism factors are unchanged, because an uncolored rooted-tree
isomorphism preserves its unique bipartite coloring relative to the root.

For different fixed roots, whole-copy pairing leaves the top pair of
edges and hence exactly Q_ij. All other doubled-tree branches sum to one.
Opposite-shore roots give Q_ij=0, as they should. No factor two is inserted
in either the Gaussian variance or the marked-return coefficient.

## 2. The isometry U and a finite explicit first nonlinear check

Use precisely the same rooted odd-degree tree shapes, with product B_e
coefficients and the same rooted automorphism normalization. The finite
one-root old frame therefore has independent standard Gaussian limits.
Removing an external edge still yields each even normalized child-Hermite
monomial exactly once. Thus U h_T=G_T is the same even-L2 to first-chaos
isometry. A two-type Gaussian space is not required when the same finite
policy is used on both shores; each shore has the identical rooted local
law. The joint two-root law retains the possibly coherent Q kernel.

For an explicit normalization test let X1_i=sum_j B_ij S_j and

    X3_i=sqrt(2) sum_j sum_(k<l; k,l!=i)
               B_ij B_jk B_jl S_j S_k S_l.

This is the three-marked-vertex tree with external-root degree one.
On bipartite full-sign support, for every C, not just a random matrix,

    Var(X3_i)=(m-1)(m-2)/m^2,
    Cov(X3_i,X3_r)=((m-2)(m-3)/m^2) Q_ir
                       for distinct same-shore i,r,
    Cov(X3_i,X3_r)=0 for opposite-shore i,r.

Also X1 and X3 have exactly orthogonal original input-spin degrees.
Deleting the output-root collision from
B[S (X1^2-1)/sqrt(2)] yields exactly the displayed X3 kernel. This checks
both the variance-one limit and the U inverse h=(G1^2-1)/sqrt(2), including
the shore restrictions which a dense-count shortcut would miss.

## 3. Tensor estimates and the two exceptional mixed contractions

Every proof using an entry upper bound replaces (N-1)^-1/2 by m^-1/2.
In the crossing-edge compression estimate, N^(I/2)m^(-E_int/2)
is bounded by a tree-dependent constant because I<=E_int and N=2m.
The fixed-root improvement is still O(m^-1/2). Collision Hilbert norms,
row absolute sums O(sqrt(N)), and bounded root-map estimates have the
same asymptotic orders. All trees and polynomial degrees remain fixed
before m tends to infinity, so tree-dependent powers of two are harmless.

The old-input/nonlinear-forest endpoint graph with no free-free parity
edge has an optional direct a--j edge and p>=2 paths a--u--j. If the
direct edge is present, this graph contains a triangle, and the bipartite
support makes the entire term zero. Otherwise a and j must be on the
same shore, and summing each intermediate label gives the same Q_aj
factor. Any additional support constraints either are consistent or make
the graph zero. Thus the surviving matrix is a constant multiple of
m^-1/2(Q^(circ p)-I), with exactly the same O(1) Frobenius bound.

In the full-return/root-hit doubled diagram, the proof that a parity edge
survives is purely combinatorial and unchanged. Its normalization changes
by only a fixed constant; the discrepancy estimate of section 1 still
makes E||J||_F^2=o(N). In the root-not-hit case the Gaussian pure-chaos
covariance-of-squares and creation inequalities require the same proper
cut bounds and old/return covariance estimates, not nonzero same-shore
matrix entries. Fourth-order Rademacher replacement and own-spin removal
likewise use those influence bounds, not irreducibility.

The marked diagonal return remains Q_ii=1. Off-diagonal terms use
sum_j|Q_ij|=O(sqrt(N)) and the O(m^-1/2) derivatives. Block diagonality of
Q only removes some of these terms. The universal odd-Schur inequality
Tr Q sqrt(sum_(k odd)w_k Q^(circ k))>=N applies to every correlation
matrix and explicitly requires no support/coherence condition.

## 4. Actual cap conversion and equal-shore spectral deletion

At fixed ||B||op<=L the unchanged marked analytic inequality gives

    liminf (1/N) sum_i E[H(X_i)|(BF(X))_i|] >= J(F,H).

For |F|+H<=1, the actual feasible means F+H sign(BF) and
-F+H sign(BF) have quadratic energy difference

    2 sqrt(m) sum_i H_i |(BF)_i|.

Independent coordinate rounding is exact because A is hollow; both
energies are in [-L(C),L(C)]. Hence

    liminf L(C)/(2m sqrt(m)) >= J(F,H).                (3)

Insert the already certified fixed Gaussian construction only after
establishing (3) for every fixed finite approximation. This yields c_*.

To remove the operator bound, restrict attention to L(C)<=K m^(3/2).
The simultaneous diagonal Grothendieck majorant satisfies
D>=A,-A and Tr D<=K_G beta(A)<=2K_G K m^(3/2). Delete all coordinates
with D_ii>T sqrt(m). At most 2K_G K m/T coordinates are removed in total.
Remove extra arbitrary coordinates from the less-deleted shore until
both retained shores have size m'. For T=2K_G K/epsilon,
m'>=(1-epsilon)m. The retained cross block is still a SQUARE full-sign
matrix, and ||A_sub/sqrt(m')||op<=T/sqrt(1-epsilon).
Principal cap monotonicity yields L(C_sub)<=L(C). Apply (3) at fixed
epsilon/L and then send epsilon to zero; the loss is (m'/m)^(3/2).
Sequences with unbounded normalized L(C) require no deletion.

## 5. Scope of the audit

Read directly: the full fresh lower-chain reconstruction; the elementary
tree moment argument; injective input Gram and collision-operator proof;
crossing-edge tensor compression; nonlinear covariance trace proof,
especially its endpoint-graph and root-hit sections; restricted marked
center update; and the complete odd-Schur proof. No new numerical
certificate or stronger Gaussian-chaos theorem is needed for this support
change. The scalar policy remains fixed before the matrix limit.

This verifies the bipartite adaptation, not an arbitrary atomless weighted
extension. General weighted powers do not reduce to signs/support masks,
and general sparse support can leave nontrivial even-edge constraints.
Neither broader extension is imported into (2).

Exact diagnostic replay:
`computations/decisive_bridge_bipartite_tree_normalization_checks_2026_09_07.py`
passes all 530 full-sign cross blocks of orders 1,2,3 and 30 sampled blocks
of order 4 (560 total). Integer sparse-polynomial Gram calculations and
independent complete spin sums check X1/X3 covariance and orthogonality.
Exact spin enumeration checks beta(A)=2L(C) and Qcap(A)=L(C). A further
2240 rational feasible-mean tests check the energy-difference identity
and that both rounded-cube energies lie in the actual cap interval.
These checks complement the diagram proof, not replace its ordered limit.
