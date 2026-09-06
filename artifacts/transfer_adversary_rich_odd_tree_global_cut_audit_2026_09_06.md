# Uniform global cuts for the selected finite odd-tree source

2026-09-06. Independent hypothesis audit for
`transfer_adversary_bounded_response_nuclear_dual_split_2026_09_06.md`.
Status: **PASS for fixed finite odd tree/forest equality diagrams.** This
is stronger than the generic random-derivative estimate with logarithms,
but its diagram hypothesis is narrower than an arbitrary centered circuit.

## 1. Precise graph property and its preservation

Represent a fixed rooted odd polynomial computation by its finite matrix
graph. Matrix transports are edges, with uniformly bounded operator norms;
fixed bounded deterministic row coefficients are vertex weights. A seed
occurrence is a mark at its corresponding graph vertex. Before label
identifications, each non-output-root vertex v satisfies

```
deg_graph(v) == number_of_original_marks(v)    (mod 2).  (1)
```

Indeed an internal odd coordinate monomial has c child arms and m own
seed marks with c+m odd. Its parent arm adds one graph edge, so
`1+c == m (mod 2)`. Unmarked transport-chain vertices have degree two;
an ordinary marked leaf has degree and mark count one. This includes
the marked cubic source S*h_2(G): two Gaussian child arms, one own mark,
and the parent arm. Hermite subtraction terms are handled termwise.

Products of branches meeting at the output root need not satisfy (1)
at that root: it is declared a boundary port. Every other vertex does.
Expanding exact Boolean seed equalities merges vertices and adds mark
multiplicities. Graph degrees, including loops counted twice and all
parallel edges, add modulo two. Thus (1) survives every coarsening.

Odd mark blocks are the free original Walsh labels; even blocks are
summed. Consequently a nonroot vertex without a free Walsh label has
even graph degree. Root/free coincidences are permitted, with separate
identity boundary ports attached to the same underlying vertex.

## 2. Every bridge-forest leaf has a boundary

Declare the output root AND every free Walsh label to be boundary ports.
Suppose a leaf component after deleting graph bridges had no boundary.
It would contain no exceptional output root and no odd-mark vertex.
By (1), all its vertices would have even graph degree. But the sum of
their graph degrees is twice the number of its internal edges plus
ONE exiting bridge, an odd number. This contradiction shows that each
bridge-forest leaf contains a declared boundary.

If the graph has no bridge, its single component already contains the
output root and free ports. The selected positive-degree source is
connected to its root; one does not discard its matrix edges and then
silently sum a detached constant component.

The already reconstructed prescribed-port matrix-graph norm theorem
therefore applies. For ANY prescribed global flattening of this boundary
tensor, its norm is bounded by the product of its edge operator norms
and vertex-weight suprema. The root can be on either side, and arbitrary
free labels can share its side. Because the number of edges, weights,
and equality patterns is fixed, the bound is a constant independent of n.

## 3. Exact source restrictions do not introduce logarithms

The finite inclusion-exclusion implementing exact source restrictions
is carried out at graph level. A coarsening that identifies two declared
distinct free labels vanishes after the final free-distinct projection.
All remaining coarsenings retain their free-label count and the parity
property above. A final same-side distinctness restriction is a coordinate
projection; a cross-cut restriction is a rectangular pinching complement,
costing at most two. There are finitely many such restrictions at fixed
degree.

Thus every selected exact positive original-degree block of a fixed
finite odd-tree/forest source, INCLUDING its single B transport, has
uniformly bounded global coefficient cuts. This verifies the uniform-
global-cut premise used for X in the dual-split theorem whenever its
actual source has this diagram representation.

The same global bound for the untransported k-fold exact noise product
follows by tensoring its bounded global factors and identifying their
common output root. Whole branches on the root side contribute their
bounded row Hilbert norms; the remaining branch cuts are global proper
cuts. Exact distinctness again changes only degree-dependent constants.
The subsequent small-entry weighted transport uses the separate
two-factor transport lemma, not a parity assertion about an arbitrary
newly transported even polynomial.

## 4. What this verification does not cover

An arbitrary centered polynomial transport may have an even node whose
own mark count does not satisfy (1). Centering alone is not this parity
hypothesis. The generic derivative theorem remains the safe polylogarithmic
bound for such a computation. Likewise, a bounded trigonometric response
has infinitely many original degrees: its exact projections have not
thereby acquired the finite selected-diagram structure above.

This audit supplies uniform global cuts for the specified finite odd-tree
source. It does not settle the rich returned-query comparison involving
the entire bounded BC, nor assert unlimited-depth feedback closure.

## 5. Two related analytic checks

The root's explicit continuous-moment counterexample
`transfer_director_continuous_moment_approximation_counterexample_2026_09_06.md`
passes independently. With W=X^3 and
`h(w)=sin(sqrt(3)|w|^(2/3)/2-pi/6)`, the k=2m moment has gamma factor
Gamma(3m+1/2) and phase
`-pi/6+(3m+1/2)pi/3=m*pi`; odd moments vanish by parity. Its nonzero
L2 norm and exact polynomial orthogonality follow. It is not a counterexample
to a direct Boolean-Stein argument or to the actual first-marked theorem.

The root's symmetrized cubic chain remainder also passes with an explicit
constant. For `h(t)=f(v0+t*d)`,

```
(h(1)-h(-1))/2 - (h'(1)+h'(-1))/2
 = -(1/4) integral_(-1)^1 (1-t^2) h'''(t) dt.
```

Hence its absolute value is at most `||D^3 f||infinity ||d||^3/3`.
Since Delta_a Z is independent of S_a, conditional averaging permits
this identity inside `E Delta_a Z Delta_a f(V)`. There is no quadratic
remainder after that averaging. High-influence coherent queries can
still leave only an O(n^-1/2) entrywise gain; this cubic identity alone
does not justify nuclear accuracy. The dual split addresses that
dimension issue by a different argument.

The seed's bounded trigonometric global-cut lemma was also read and
reconstructed in full. The entrywise remainder of exp(itJ)-1-itJ has
absolute row and column sums at most `(t^2/2)||J||op^2`, proving the
dimension-free quadratic matrix bound. Its exact finite-difference
cover formula follows by the alternating subset sum. The embedding of
the full seed tuple into its covered subtuples is injective precisely
because they cover every differentiated label, so the tensor-product
operator argument is valid. Fourier projection then identifies the
claimed global coefficient cuts. This supplies polylogarithmic cuts,
not the uniform selected-odd-tree cuts proved above and not small proper
cuts of an arbitrary returned bounded response.

Read and reran
`computations/transfer_seed_trigonometric_cut_checks_2026_09_06.py`:
192 entrywise exponential matrix cases, 3200 exact Boolean cover scalar
cases, and 96 modified derivative matrix cases passed. These checks
test the stated finite algebra, not rich feedback closure.
