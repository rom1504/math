# Independent audit: bounded-amplitude weighted marked lower bound

2026-09-07. **PASS**, with the unchanged scalar lower certificate as an
explicit dependency. This audit independently reconstructs the sensitive
weighted changes; it is not a fresh numerical replay of that certificate.

Target source: `principle_invent_2026_09_07_weighted_tree_extension.md`.
The old input Gram and full nonlinear return sources were reread through
their actual combinatorial proofs, not accepted from prior verdicts.

## 1. Off-diagonal old-input Gram

Two input trees with equal marked vertex count q have 2(q-1) edge
occurrences. Every edge occurs once or twice. Their parity graph therefore
has an even number of edges; off the output diagonal it is nonempty.
The apparent old exceptional graph consisting only of the root edge is
impossible. Doubling the Frobenius expression leaves a single edge with
at least one unglued free endpoint. All other weighted factors split into
bounded unary functions after the other labels are fixed. Thus the actual
weighted bilinear norm, not a fictitious sign-only norm, supplies the same
O(N^(-1/2)) density gain. Empty-parity diagonal patterns are doubled trees;
their leaves sum by row variance. Output-root collision estimates retain
their exact label count and bounded-amplitude coefficients.

## 2. Exceptional full contraction

I reconstructed the common-tree lemma separately before reading its final
write-up. With roots a,j and no free-free parity edge, the optional a-j
edge must occur by odd edge-count parity. Since j is absent from the old
input tree, all j edges belong exclusively to the new forest tree.
Every a-u path edge must belong to the old tree, or it closes a triangle
in the new one. The common forest has p+1 components; a and the p path
middles must occupy distinct components, since a common path would close
a cycle in one of the two trees. Hence each component has exactly one
anchor. This is a tree theorem, not an assumption on a variance profile.

Summing its doubled edges gives bounded POSITIVE unary anchor weights.
The remaining matrix is diag(f_0) times

    B circ (B diag(f_1) B) circ ... circ (B diag(f_p) B), p>=2.

Each inner matrix is PSD with bounded operator norm and bounded row
Euclidean norms. Two factors and row Cauchy--Schwarz give absolute row
sum O(N^(-1/2)), and Frobenius O(1). This retains arbitrary squared-edge
weights and is stronger than replacing them by approximately constant
scalars. The root-hit doubled-diagram argument separately keeps a single
free parity edge and gives E||J||F^2=o(N).

## 3. Remaining modules and normalization

The existing global/fixed-root tree cuts depend on bounded operator norm,
row norms, and maxentry O(N^(-1/2)), so changing the entry bound from 1 to
fixed K changes their fixed constants only. Whole-branch matches use the
actual old-input covariance just checked. Partial matches, collision
deletion, bounded-degree hybrid replacement and smoothing retain their
finite-before-order limits.

For uniform row error eta_N=o(1), normalize Q=B^2 by its diagonal to a
correlation matrix Qhat, not by falsely rescaling B symmetrically to exact
row norm. At fixed L, ||Q-Qhat||op=O_L(eta_N); fixed Schur polynomials
change by o_op(1). The trace/square-root variance inequality therefore
changes by o(N). This validates the final scalar certificate passage.

Finally the director's variance-completion lemma removes the operator cap:
delete epsilon N heavy coordinates, rescale mildly, complete the variance
with a small independent symmetric perturbation, apply the fixed-operator
bound, then send epsilon to zero. The invention agent independently checked
the exact deficit realization and all error normalizations.

## Verified consequence and strict scope

Every symmetric hollow real W_N with a FIXED entry bound and uniform
row squared norms N-1+o(N) satisfies

    liminf width(W_N)/N^(3/2) >= 0.4333221116640807.

Width is half the full energy range; the same bound follows for cap.
No mixing, finite species, spectral flatness or sign pattern assumption
is present. The fixed scalar certificate and original analytic marked
policy remain dependencies. This extends an existing lower theorem; it
does not improve its constant.

In particular, an amplitude profile cannot beat the certified constant
merely by remaining bounded and row-regular. Equality of its optimized
value with the full-sign optimum above that constant is NOT proved.
This is not sign rounding, a recurrence, or convergence. Allowing entries
to grow as sqrt(N) admits the explicit five-cycle-block .4 counterexample;
no claim covers that case.
