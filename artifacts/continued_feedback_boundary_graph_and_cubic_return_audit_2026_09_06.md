# Fresh boundary-graph audit and exact cubic coherent-return subtraction

Date: 2026-09-06. Independent feedback-agent reconstruction of
`continued_director_boundary_graph_feedback_cuts_2026_09_06.md`.
The main high-degree estimate passes. Section 4 gives a stronger cubic
classification; the audit agent has independently read and passed its
full proof, including derivative multiplicities and exact exclusions.

## 1. Primary-source and prescribed-boundary audit

I read Theorem 11 and Lemma 15, including their complete proofs, in
[Mingo--Speicher, Sharp Bounds for Sums Associated to Graphs of Matrices](https://arxiv.org/html/0909.4277).
Theorem 11 bounds the input-output graph operator by the product of edge
norms; it does not bound arbitrary partial transposes automatically.
Lemma 15 supplies an acyclic input-output modification through identity
vertex splits. The director's extra boundary argument is needed.

That argument passes: adjoining the two hubs and the edge between them
removes every bridge under boundary-leaf coverage. Deleting either hub
leaves a connected graph. This remains true after internal identity
splits. Hence an unused edge at a hub can be extended to an ear returning
elsewhere, by a path toward the other hub avoiding the first. Neither
hub need split. A degree-two port enters an ear with both its edges, so
it never subsequently splits. Removing the hubs leaves exactly the
prescribed boundary orientation and no internal source or sink.

Thus I find no hidden dimension factor from a closed cycle or from the
chosen two-sided boundary partition. Vertex weights and final rectangular
block pinchings have the stated operator norms.

## 2. Exact equality bookkeeping and the three-mark bridge bound

The equality partition must include the distinguished root label `j`,
even when no original own-spin factor marks that root. This is necessary
for diagonal `Q` contributions. It is a clarification of the exact
bookkeeping, not an additional restriction on the theorem.

Every original non-root vertex has degree parity equal to its seed-mark
parity. Identification preserves this fact, with loops contributing two
to degree. The root's parity is allowed to differ. After even seed
blocks are summed, every internal non-root vertex has even degree.
Identifying an internal block with a free block leaves its parity odd;
identifying two internal blocks leaves it even. Thus internal/free
inclusion-exclusion preserves the required parity description. Distinct
free labels are excluded only by the final coordinate projections and
cross-cut block pinchings.
If an inclusion-exclusion term transitively merges two prescribed free
vertices, the final free-distinctness projection kills it; it is not
retained as a graph with the original number of separate boundary vertices.

For a bridge side away from the source root, every primitive center
inside supplies its own direct root edge across the cut. There can
therefore be at most one center. If there is one, any seed vertex of
another primitive inside would supply a second crossing edge. The side
contains only the selected primitive's three seed marks. With no center,
it contains only one leaf. This remains true with identified vertices
because parallel edges are retained.

An uncovered non-root leaf component is impossible by even-degree
parity and its unique exiting bridge. If the root component is an
uncovered leaf, all free marks lie beyond its one bridge and hence
there are at most three. This proves the director's `q>3` coverage
criterion without discarding any exact Walsh corrections.

The exact finite combinatorial audit enumerates every vertex partition
of the `Y^3` diagram (115975 partitions, including root identifications),
and all partitions of `Y^2G`, `YG^2`, `G^3`, `SY^2`, `S^2Y`, `S^3`,
`Y^2`, and `SGY`. It tests degree parity, every away-bridge mark count,
boundary coverage, and the cubic exception below. The enumeration is
adversarial evidence; the preceding argument handles arbitrary fixed
polynomial degree.

## 3. Two immediate sharpenings

First, the exact returned degree-two tensors also have all proper cuts
`O_{P,L}(epsilon)`. Any source monomial contributing to an even Walsh
degree contains an even number of odd primitives. The root degree
parity then agrees with its seed-mark parity as well. Consequently ALL
internal vertices are even, including the root; a boundary-free leaf
component anywhere is impossible. This covers degree two directly.

Second, declare the source root itself an additional boundary port.
Every leaf component is now covered for every positive Walsh degree.
The same graph theorem gives an `O_{P,L}(1)` bound on every GLOBAL
proper cut of the exact source tensor, whose slots are the root and
its positive-degree Walsh labels. If a root label coincides with a free
label, equality embeddings or block pinchings handle the repeated
boundary port. This strengthens the earlier polylogarithmic global-cut
module for this particular canonical history.

The cubic exception is real. On a conference sequence with `Q=I`,
take `P(W)=Y`. Then `BP=BY=D`. In a symmetric ordered degree-three
coefficient convention,

`K_i(i,a,b)=(sqrt(2)/6) B_(ia) B_(ib)` for distinct `i,a,b`.

The row of its one-slot-versus-two-slot flattening indexed by `i` has
squared norm `(1-sum_a B_(ia)^4)/18`. This stays bounded away from
zero as `epsilon->0`. Thus omitting the cubic qualification would be
false on actual hollow signings.

## 4. Exact cubic subtraction theorem

Retain the director's hypotheses: `B` is symmetric and hollow, its
squared row norms are one, `||B||<=L`, and `max|B_ij|<=epsilon`.
Write

`W=(S,G,Y,QS,QD)`, `G=BS`, `D=S h_2(G)`, `Y=BD`.

For a fixed polynomial `P` in these five FORMAL coordinates define
the deterministic vectors

`u_j=E[partial_3 P(W_j)]`,
`v_j=E[partial_5 P(W_j)]`.                                  (3)

Let `R_i` be the exact degree-three Walsh component of

`B [ P(W) - D_u Y - D_v QD ]`.                              (4)

Then every proper flattening of its symmetric degree-three coefficient
tensor satisfies

`max_i ||K_(R_i)||_(cut,op) <= C_(P,L) epsilon`.              (5)

The derivatives are with respect to the displayed formal coordinates;
one does not first simplify accidental algebraic relations among W.
This convention permits a nonunique coherent/noise decomposition, but
the estimate holds for every fixed polynomial representation.

### Proof

Expand one source monomial into labeled primitive diagrams, retaining
all multiplicities. An uncovered degree-three equality diagram has
exactly one uncovered leaf component: the root component. It contains
no free mark. The other side of its sole bridge has three free marks.
By Section 2 it must consist of one distinguished primitive D-star,
with its center and its two leaves all distinct and free. Its direct
root edge is `B` if the selected factor was Y and `Q` if it was QD.
No other primitive can meet any of these three vertices: such a meeting
would supply a second edge across the alleged bridge.

All other factors are therefore a closed root-side expectation diagram.
Conversely every nonzero expectation diagram of these remaining even
factors has even degree at every vertex, including the root. It is
connected and has no bridge. Attaching a free distinguished D-star to
it therefore produces precisely the uncovered configuration in question.
Choosing the distinguished Y factor in all possible ways gives exactly
the formal derivative `partial_3 P`; choosing QD gives `partial_5 P`.
The three-free-mark component is unique, so there is no extra count
from selecting two different free stars in the same uncovered diagram.

Initially restrict this root-side expectation's internal labels to be
different from the three displayed free labels. Removing these
restrictions changes only diagrams in which a free star vertex is
identified with a root-side vertex. Every such correction is covered:
an identification with the root component puts a free mark there;
an identification creating another connection destroys the sole-bridge
separation. The non-root leaf parity argument handles every other leaf.
The same statement covers equality of the star center and the source
root when the direct matrix is Q.

After those restrictions are removed, the closed expectation factors
exactly as (3) times Y or QD. Hence all uncovered terms in P cancel
against (4), term by term with their correct source multiplicities.
What remains is a finite sum of covered diagrams. Apply the prescribed
boundary graph theorem to each, with the outer weight `B_ij`, obtaining
an epsilon factor. Final exclusion of coincident free labels changes
the bound only by a fixed number of norm-at-most-two pinch complements.
This proves (5).

The coefficients u and v are uniformly bounded. For completeness,
`Cov(D)` has diagonal `1-sum_a B_(ja)^4` and off-diagonal
`2 B_(ij)^2 Q_(ij)`. The Schur multiplier of the correlation matrix Q
is operator-contractive, and the symmetric nonnegative matrix `B^circ2`
has row sums one. Thus `||Cov(D)||<=3`. Every coordinate of W has
bounded variance and Boolean degree at most three. Fixed-degree Boolean
hypercontractivity gives uniform fixed moments and hence
`||u||_infinity+||v||_infinity<=C_(P,L)`.

## 5. Precise consequence and limits

The exact cubic return is now identified modulo a small-cut tensor:

`P_3[B P(W)] = B D_u B D + B D_v Q D + R_3`,

where `R_3` obeys (5). The degree-one part can simply be retained as
its exact Boolean Walsh matrix applied to S; its operator norm is
bounded by the source global-cut result and `||B||`.

This is a finite-polynomial structural theorem. Small proper cuts
do NOT mean small L2 norm, nor independence from every other cubic
field: equal-degree covariances must still be retained. A bounded
response requires the already ordered approximation passage; a positive
innovation variance or an improved original constant is not implied.

Code and results:

- `computations/continued_feedback_boundary_graph_partition_audit_2026_09_06.py`
- `computations/results/continued_feedback_boundary_graph_partition_audit_2026_09_06.json`

## 6. Fixed finite coherent-family corollary

The same proof applies to a fixed finite coordinate family

`W=(M_1 S,...,M_a S, N_1 D,...,N_b D)`,

where every deterministic matrix has bounded operator norm. A linear
primitive has one direct root edge and one marked seed; a cubic primitive
has one direct root edge N_beta and the same two B edges and three marks
as before. Therefore the bridge count, global-cut estimate, and cubic
classification are unchanged. In particular, for fixed polynomial P,

`P_3[B P(W)] = sum_beta B D_(u_beta) N_beta D + R_3`,

`(u_beta)_j=E[partial_(a+beta)P(W_j)]`,

with all proper cuts of R_3 bounded by `C epsilon`. All higher Walsh
degrees have the same small-cut bound, and the exact first Walsh matrix
can be kept literally. The matrices `B D_(u_beta) N_beta` remain bounded
in operator norm, so this is a finite coherent cubic-state closure.

This corollary does not remove the new high-chaos fields from a subsequent
nonlinear response. It only classifies its coherent degree-three return
when the response's displayed coordinates belong to this fixed family.

## 7. Fresh check of the next-return open-mark assembly

I also read both complete audit notes
`continued_audit_next_return_open_mark_hall_2026_09_06.md` and
`continued_audit_boundary_graph_and_next_return_2026_09_06.md`.
Their fixed-polynomial full-contraction conclusion is supported by the
following explicit version of the delicate unsplit factorization.

Let the left source be `(A-EA) Z`, at a fixed positive even Walsh degree
of A and a fixed odd degree of Z. Let the right source be a polynomial
in old G,Y, and keep every unmatched right mark. In an unsplit diagram,
put all right primitives meeting Z into the first group. Any right
internal label shared with the other group must be retained as an OPEN
index K, not summed in either group separately.

After fully contracting Z, the first group gives a tensor
`S_(a,b,K,U_S)`. The small proper cut / fixed-slot / exceptional full
degree-three comparison yields

`sup_(a,b) ||S_(a,b,.,.)||_F <= C epsilon`.

This is a Frobenius bound INCLUDING K, stronger than an estimate at each
fixed value of K. The remaining group gives
`Gamma_(a,b,K,U_Gamma)` by fully contracting the positive-degree A
tensor. Its selected right diagram has a bounded global cut from A's
contracted labels to all remaining indices `(b,K,U_Gamma)`. Composing
that map with A's root tensor therefore gives

`||Gamma||_F <= C sqrt(n)`.

The full diagram is a contraction over K of S and Gamma, with repeated
remaining labels retained by equality projections. Pointwise Frobenius
Cauchy--Schwarz followed by the two-root sum now gives

`||H||_F^2 <= C epsilon^2 ||Gamma||_F^2 = O(1)`.

Thus there is no uncharged factor from a sum over labels that had merely
been pinned. In the split case, the old Y is partially merged with A
and another old primitive with Z; these vertex-disjoint merges supply
two small gains while all shared hyperedge labels remain explicit.

Finally apply B and `B D_a B` on the two root axes and only then
restrict their output roots to the diagonal. These are bounded operators
followed by a norm-one projection. The full-contraction squared norm
divided by n vanishes. Local collision errors are handled in the same
order using their positive-degree root Frobenius error and the selected
right-diagram global cut. The Bc0 part is more direct: its higher-degree
right tensor has a small proper cut by the director's theorem.

This verifies the fixed-polynomial next-return criterion, with equal-
degree covariance retained. It does not skip the remaining ordered
bounded-response/variance passage or assert positive innovation.
