# Independent audits: typed rank-two certificate and regular-support weave

2026-09-07. Both displayed mechanisms pass first-principles reconstruction.
Their remaining scopes differ: the regular-support extension gives a new
actual low-bilinear-cap bridge, whereas the rank-two typed criterion still
needs a growing-order inequality verified before yielding a new upper.

## Rank-two typed certificate: PASS

Audited source:
`artifacts/flatify_director_rank_two_typed_certificate_2026_09_07.md`.

For each fixed physical spin row, the random signed column permutation is
exactly a uniformly permuted magnitude word with independent fresh signs.
Replacing that word by iid samples conditioned on its exact type is an
identity, not an approximation. On the conditioning event the inserted
potential sum is k E_rho phi. Dropping the indicator from a nonnegative
expectation consequently gives precisely the negative potential term and
the factor 1/p_k(rho) in the source.

There are k-2 coordinates used in cross pairs and two unused loop
coordinates per row. The latter give exactly (E exp(phi))^2. A useful
independent consistency check is potential gauge invariance: adding a
constant c to one row potential contributes (k-2)c from cross factors,
2c from loops and -kc from the conditioning correction. These cancel.
Omitting the loop factors would fail even this elementary test.

The threshold is correctly 2tqN: Q(C_cross)>=qN^(3/2) means
|x^T C_cross x/sqrt(N)|>=2qN. Both signs have the same distribution by a
fixed bijection of column signs: for every unordered fibre edge flip both
coordinates at one endpoint. Different pairs use disjoint coordinates.
The physical preimage counts c_k include all Boolean rows; any duplicate
global-sign representations merely make the valid union bound looser.
Potential infima may depend on the entire type tuple because the argument
is a pointwise upper bound before a finite sum. Nonattainment is handled
by the strict certificate margin.

The type overhead estimates also pass. A type with positive integer
counts n_a is a mode of the multinomial with rho_a=n_a/k: the probability
ratio of any other count vector to this one is at most one, factor by
factor. There are at most (k+1)^s count vectors on its s-symbol support,
so p_k(rho)>=(k+1)^(-s). Distinct positive integer magnitudes have squares
summing at least 1^2+...+D^2, giving D=O(k^(2/3)). The square-partition
generating function bound in the source independently yields
log |R_k|=O(k^(2/3)). Multiplying the logarithmic errors by m=k/2 still
gives o(N), uniformly over all realizable rows.

Filling m blocks of order k by O(k^(3/2))-cap signs costs
O(m k^(3/2))=O(N^(5/4)). These are genuine full signs. No independence of
their cap events from the cross construction is required: deterministic
choice followed by triangle inequality suffices.

No gap found. The unproved item remains precisely the large-k maximum
over realizable type tuples in the certificate, not the finite reduction.
Independent column signs erase an outer scalar signing, so even a successful
certificate is not automatically a selected-actual-child transfer theorem.

## Regular-support direct-E extension: PASS after one direction correction

Audited source:
`artifacts/flatify_independent_2026_09_07_regular_support_weave.md`, with the
fixed-depth row theorem reconstructed in
`artifacts/decisive_audit_standalone_direct_E_upper_2026_09_07.md`.

For d-regular loopless support on m fibres and k retained rows, each
undirected edge is counted twice, giving exactly

    D_sigma=2(mdk-sigma x^T W x).

Averaging its sign in exp(-t D_sigma/(2k)) therefore produces the stated
folded kernel, without another factor of two. There are no self-loop
coordinates or deleted-coordinate losses in this support version.

The graph tensor contraction is valid for every finite loopless graph.
One elementary proof avoids any ambiguity about a generalized Holder
statement: merge two adjacent tensors, contracting ALL their common
indices, and use Frobenius submultiplicativity for the resulting matrix
product. Continue until each connected component is a scalar. Its modulus
is at most the product of all original Hilbert norms. Contracting all
common indices prevents internal loops during merging. Finite Gaussian
feature truncations give the result first in finite dimension; the
nonnegative even-feature expansion permits passage to the full kernel.
The squared norm of a permutation-averaged row tensor is exactly the
signed-permutation orbit quantity P_t, hence its norm is L_t as required.

The row theorem's o(d) error is uniform in terminal Hadamards and does not
depend on m. Multiplication over independent fibres yields o(md), including
m=2d. The exponent is t gamma+p log2+E, so gamma must approach
-(p log2+E)/t from BELOW for a negative exponent. The initial draft said
"above"; this was reported to its author. Approaching from below gives the
same limiting cap constant C from above and changes no endpoint.

For K_(d,d), the matrix is [[0,B],[B^T,0]], with a full square sign bridge
of order n=dk. Its quadratic cap equals the bilinear cap of B exactly.
The general bound C sqrt(d/m) N^(3/2), at m=2d and N=2n, becomes
2C n^(3/2), not C or sqrt(2)C. Thus the audited rational point gives

    limsup R_n/n^(3/2) <= 2U < .987216188.

At fixed recursive depth the previously established multiplicatively
dense Hadamard orders give multiplicatively dense n=d floor(pd).
Restricting both row and column sets of B cannot increase its bilinear
cap: average independent unbiased omitted spins. This supplies every order
without a growing-depth premise.

This is a genuine correlated full-sign bridge operation. The inequality
Q(parent)<=Q(A)+Q(D)+beta(B) is still too weak to pay actual-child shells,
so no recurrence or convergence is inferred.
