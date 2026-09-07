# Querywise placement does not commute with the actual maximum

2026-09-07. **Proved scoped falsifier.** This is an exact quantitative
reformulation of the archived reciprocal-weave matching floor, not a new
impossibility theorem for arbitrary global edge changes.

Let m>=3 admit a full Hadamard F of order m, and let S be any symmetric
full-sign macro matrix. Use m physical fibres of size m. At each fibre i
one may choose a permutation pi_i of the m Hadamard columns. For distinct
fibres, use the actual full-sign outer-product tile

    C_ij = S_ij f_(pi_i(j)) f_(pi_j(i))^T.

Use fixed hollow full-sign internal blocks with cap O(m^(3/2)). The
resulting actual parent has order N=m^2. Consider only the query family

    X = { x_i=epsilon_i f_(a_i): a_i in [m], epsilon_i in {+1,-1} }.

Its cardinality is at most (2m)^m, with logarithm O(m log m)=o(N).

For every prescribed x in this family, there exists a choice of the
column permutations for which ALL cross energy is exactly zero. Namely
choose a directed m-cycle and send the selected column at fibre i to its
successor port. A cross tile can contribute only if the two directed
choices are reciprocal, and the cycle has no reciprocal pair. Consequently

    max_(x in X) min_pi |H_(C(pi))(x)| = O(m^(5/2)).

On the other hand, for every fixed choice of all pi_i, match the fibres
in floor(m/2) pairs and choose each matched fibre's query column to be its
partner port. Every matched cross term has magnitude m^2, and all other
cross terms involving these fibres vanish. The signs epsilon_i can make
all these contributions have a common polarity. For even m, the total is
exactly m^3/2 in magnitude before internal filling. For odd m, leave one
fibre unused in the matching; choose its column arbitrarily and average
its global sign, together with independent common signs on the matched
components, to kill unpaired contributions.

The internal energy is invariant under each fibre's global sign. Choose
the common matched-edge polarity to agree with the sign of the total
internal energy, and then average the component signs. This gives

    min_pi max_(x in X) |H_(C(pi))(x)|
       >= floor(m/2) m^2
       = (1/2-o(1)) N^(3/2).

Thus the actual simultaneous maximum has a leading-order gap from the
query-by-query favorable placement, even though the query family has
entropy o(N), much smaller than the full Boolean cube. A small query
count at this scale is not a sufficient reason to exchange the two
optimizations.

The matching mechanism and its half floor are already archived. The
point recorded here is the precise min--max and query-cardinality boundary
for attempts to reuse the new query-dependent port-defect compiler. It
does not say that a sub-half construction with different tiles is
impossible, and it does not restrict the user's globally rewritten-parent
scope or selectable optimal-child choices.
