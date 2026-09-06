# Independent audit: boundary graph cuts and the next finite return

Date: 2026-09-06. Full-read reconstruction of
`continued_director_boundary_graph_feedback_cuts_2026_09_06.md`, with
the fixed-polynomial next-return consequence below. The theorem does
not supply a positive innovation bound or a new universal constant.

## 1. The prescribed-side graph lemma passes

For a connected matrix graph whose every bridge-forest leaf has a
boundary, each nontrivial flattening of its boundary tensor has norm
at most the product of the edge operator norms. A bounded weight at
any vertex multiplies this bound by its supremum norm.

The multi-input/multi-output orientation is essential: a partial
transpose of an already bounded operator would not prove this claim.
Attach separate identity ports for all prescribed boundary indices,
then temporary hubs s,t and the edge st. Every original bridge has a
boundary on both sides, so the augmentation has no bridge. A hub and
its incident auxiliary edges are used only to choose orientations;
they are not summed matrix indices and do not identify boundary labels.

The standard identity-edge ear modification can be made without
splitting s or t. Their deletions leave connected graphs, so an unused
edge at either hub has a return path to the current graph avoiding
that hub. No cycle ear based there is needed. The degree-two boundary
ports need no splitting either. Once the auxiliary edges are deleted,
the resulting acyclic network has exactly the prescribed boundary
directions and no new internal source or sink. This repairs the one
detail not supplied by an arbitrary input-output modification.

The matrix norm bound follows by composing edge maps and equality
partial isometries in topological order. Weighting one equality map
gives precisely the supremum norm of the weight. The classical norm
and ear ingredients are [Mingo--Speicher, Theorem 11 and Lemma 15](https://arxiv.org/abs/0909.4277),
which were independently retrieved and read. The prescribed-port
extension and hub-preservation argument above were reconstructed here.

Multiple identity boundary ports may attach to the same original
vertex. They enforce equality of those exposed indices by the same
norm-one gates. This version is useful when a root index also occurs
as a free seed index.

## 2. Exact marked diagrams and the parity test

For `W=(S,G,Y,QS,QD)`, the primitives are exactly the root mark, a
B/Q root-to-seed edge, or a B/Q root-to-center edge followed by two B
edges and the three distinct marks at center and leaves. The degree-
three formula is exact because unit row norm cancels the h2 constant
and hollowness excludes the marked seed from its two Gaussian slots.

Products on the Boolean cube are expanded by equality partitions of
their original seed marks. Odd blocks are free Walsh labels; even
blocks are summed. At every nonroot vertex graph degree has the same
parity as mark count. This remains true when centers and leaves merge.
Multiple edges must be retained for this parity and bridge argument.

A nonroot bridge-forest leaf without a free mark would have all even
degrees but one exiting edge, a contradiction. If the root leaf lacks
a free mark, all free marks lie beyond its one bridge. Every primitive
center has its own direct root arm. A component beyond a single bridge
therefore contains at most one primitive's three original mark
occurrences (or just one linear leaf). It cannot contain more than
three distinct free labels. Hence degree q>3 forces boundary coverage
of the root leaf as well.

An exact inclusion-exclusion qualification is necessary: simultaneous
internal/free identifications can transitively merge two declared free
vertices. Such terms vanish under the final free-distinct projection
and are discarded. All remaining coarsenings contain at most one free
vertex in each equality class; they preserve q and the parity argument.
The final exclusions between free labels are coordinate projections
or rectangular pinching complements, with degree-dependent norms.

Thus every fixed polynomial P satisfies

```math
\max_i\|[B P(W)]_{i,q}\|_{\mathrm{proper\ cut},op}
      \le C_{P,L}\max_{i,j}|B_{ij}|,\qquad q>3.          (1)
```

This controls individual exact equality diagrams, not only their sum,
and makes no false small-cut assertion about QD itself.

## 3. A useful stronger global-cut corollary

Declare the source root itself an additional boundary port. Then the
root bridge-leaf is covered automatically, at EVERY positive Walsh
degree. The nonroot leaves are covered by the same parity argument.
The prescribed-port graph lemma therefore bounds every global cut of
each individual coherent equality-diagram kernel by a degree-dependent
constant. Root/free coincidences use two ports on the same vertex.

The positive-degree qualifier is essential: at degree zero there need
not be two nonempty boundary sides, and a constant root vector may
have Euclidean norm sqrt(n). This corollary supplies bounds for the
selected/coarsened coherent diagrams that could not previously be
inferred from bounds on whole polynomial sums.

## 4. Completion of the fixed-stage higher-degree full contractions

Let

```math
C=c_0+A(W)Z+R_{\ge2},\quad a=EA,
\quad \eta=B[(A-a)Z+R_{\ge2}],
\quad L=Bc_0+BD_aB r(G,Y).
```

All expressions here are at a fixed polynomial stage after the already
audited local noise-collision surgery. The surviving eta degrees are
at least five. For a degree-p eta tensor K and a degree-q L tensor J
with q>p, the desired criterion is

```math
\frac1n\sum_i\|K_{i,p}\star_p J_{i,q}\|_F^2=o(1).       (2)
```

For Bc0 this follows immediately from (1): the full p-contraction is
a proper cut of the RIGHT degree-q tensor. Its square is bounded by
C/n times the left tensor's squared norm, whose root average is bounded.

For the return BDaB r, the at-least-two-noise part was proved in
`continued_audit_next_return_open_mark_hall_2026_09_06.md`: open-mark
Hall matching gives a source two-root Frobenius norm O(polylog n),
and B tensor BDaB followed by diagonal-root restriction preserves it.

For the remaining one-noise term X=(A-a)Z, expand the right source r
into old G/Y primitives, retaining its extra free marks U. If an old
Y splits between the centered coherent factor and the noise, merging
(A-a,Y) gives one small proper-cut gain. The noise has another old
neighbor; merging that distinct pair gives a second gain. Open marks
and remaining hyperedges are retained through the Hilbert-contractive
merges. This gives source row norms O(1/n).

In an unsplit diagram, collect all noise-neighbor primitives into an
open source S. Its conditional Frobenius norm is O(n^-1/2), uniformly
in the two roots and any pinned free labels: use a proper cut of the
noise, except for a completely contracted degree-three h3(G)/Y pair,
whose covariance was already explicitly bounded. The centered even
factor A-a has positive Walsh degree. Its contraction against the
remaining coherent diagram gives an open bridge Gamma with global
Frobenius norm O(sqrt(n)), by its root map and Section 3's global cut
of the right diagram. Thus S times Gamma has global Frobenius norm O(1).

Exclusions pinning a noise-summed label to a right free mark give the
same small fixed-slot influence. Exclusions identifying a noise label
with a left coherent label produce the local collision kernels already
bounded by O(n^-1/2) in row L2. Their degree remains POSITIVE and odd.
Their root Frobenius norm is therefore O(1), and the coarsened right
diagram has a bounded global cut by Section 3. They contribute O(1)
as well. Coarsenings joining two remaining free marks vanish after the
final free-distinct projection.

Consequently each full source tensor H_(a,b,U) has global Frobenius
norm O(polylog n). Apply B and BDaB on its two root axes and restrict
the output roots to the diagonal. This is a norm-one projection after
bounded operators, so its squared norm divided by n tends to zero.
Together with the quantitative local-surgery error bounds, this proves
(2) for the full eta and both pieces of the ACTUAL L at fixed degree.

## 5. What is, and is not, licensed next

The exact Boolean stable-noise criterion now applies at each fixed
polynomial stage. Proper eta contractions were already small; Section
4 supplies the full-noise contractions into higher degrees of L.
Equal-degree scalar covariances must be retained, not set to zero.
Regression against retained homogeneous eta components can be performed
after lower variance and upper row-variance cutoffs; the literal
regression residual is retained as a coherent variable.

This does not prove a pointwise variance floor, a nonzero innovation,
or a closed scalar state evolution. Components of eta with very small
variance may be dropped from the NEXT FIELD in averaged L2, but their
possibly large regression coefficients must not be treated as uniformly
bounded. For general bounded responses, polynomial approximation,
the n limit, degree limits, and variance cutoffs still need an explicitly
ordered passage. No unlimited-depth or bounded-response second-query
theorem is asserted by this fixed-stage audit alone.
