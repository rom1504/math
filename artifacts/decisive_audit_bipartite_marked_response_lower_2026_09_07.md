# The marked-response bound for square Gale--Berlekamp matrices

Status: proved at the exact bipartite full-sign scope below. The proof
reconstructs the changed diagram weights and every operator/return input;
it does not assume a general weighted-array extension.

For C in {+-1}^{m x m}, put

    beta(C)=max_{u,v in {+-1}^m} |u^T C v|,
    G_m=min_C beta(C),
    B=(1/sqrt(m)) [[0,C],[C^T,0]],      N=2m.

Then B is symmetric and hollow, every row has squared norm one, and
max_ij |B_ij|=m^(-1/2)=sqrt(2/N). The conclusion is

    liminf_m G_m/m^(3/2)
      >= 2 c_*
       = .866644223328161506831625857795158693117269296067386826213992
       > .8666442233281615,                                  (1)

where c_* is the certified rational lower enclosure in
`decisive_audit_certified_minimum_width_lower_2026_09_07.md`.

No assertion of a Gale--Berlekamp limit, or equality with a quadratic
signing optimum, follows. The coefficient is a lower bound for real
Boolean row/column switching, not for a continuous-phase version.

## 1. The exact normalization and finite feasible comparison

For z=(u,v), the hollow quadratic energy of B is

    H_B(z)=z^T B z/2=u^T C v/sqrt(m).

Its interval is symmetric, since replacing u by -u negates the energy.
Thus

    Q(B)=W(B)=beta(C)/sqrt(m).                               (2)

For any random vectors F,H with H>=0 and |F|+H<=1 coordinatewise,
the literal cube vectors F+H sign(BF) and -F+H sign(BF) have energy
difference

    2 sum_i H_i |(BF)_i|.

Hollowness licenses independent rounding. Equation (2) gives the exact
finite inequality

    beta(C)/(2m^(3/2))
      >= (1/N) sum_i E H_i |(BF)_i|.                         (3)

It remains to establish the same marked-response lower limit for the
right side of (3). That is done first for fixed ||B||op<=L.

## 2. What replaces deletion of squared signs

Fix any finite old tree/forest construction from the lower theorem.
Use the same normalized coefficient products, but use the actual B
entries: each edge contributes C_ab/sqrt(m) on opposite shores and zero
on the same shore. Equivalently, keep unnormalized edge factors C_ab
and normalize each edge by sqrt(m).

Every relevant connected overlay of these trees has the following exact
property. If the overlay is not bipartite it contributes zero. Otherwise,
fixing the shore of one output vertex forces the shore of every other
vertex in that connected component. Each free label then ranges over
exactly m coordinates, subject only to its fixed shore and finite
injectivity exclusions. When several components occur, each contains a
fixed output vertex; their color choices are likewise fixed.

On such a consistent overlay, paired C factors really do cancel:
C_ab^2=1 on all its allowed labelings. A doubled edge still enforces the
shore coloring, so it must not be silently dropped before checking
consistency. After that check, it has no further weight.

For a component with v free vertices and finitely many forbidden labels,
the allowed unconstrained count is m^v and the injectivity correction is
O(m^(v-1)). Fixed opposite-shore coordinates are automatically distinct.
Constants throughout depend only on the fixed construction and L.

This is stronger than a generic bounded-amplitude assertion: for arbitrary
weighted arrays, squared-edge weights can remain nontrivial functions of
several free labels. In particular, a no-free-free-edge diagram need not
reduce to an entrywise power of B^2. That broader extension is not used.

## 3. Old Gaussian frame: leading moments and marked Gram

For an external-root tree T with d marked vertices and rooted automorphism
count a(T), use the field

    X_(T,i)=a(T)^(-1/2) sum_injective product_edges B * product_marks S.

Its external root has degree one and every marked vertex has odd degree.
In a fixed moment with D marked occurrences, a nonzero sign expectation
has at most D/2 free labels. Fewer labels lose O(1/m), since the edge
normalization is m^(-D/2) and each fixed shore contains m labels.

In a leading pattern every free label occurs exactly twice. An edge has
multiplicity at most two: any occurrence uses one occurrence of a free
endpoint, which has only two occurrences in the entire pattern. Thus a
surviving parity edge has multiplicity exactly one, not a hidden higher
odd power of C.

The parity graph is Eulerian. If nonempty, it has a free-free edge.
Condition on all other labels. If the forced colors are inconsistent the
term is zero; otherwise its remaining two-variable factor is one C_ab
times separate bounded unary functions. Hence its normalized sum is at
most beta(C)/m^2. Fixed ||B||op<=L implies beta(C)<=L m^(3/2), so this is
O_L(m^(-1/2)). Injectivity is imposed by unary exclusions, with the
remaining collision correction O(1/m).

If the parity graph is empty, the existing complete-copy-pairing argument
uses only connectivity, the exact vertex/edge count, and the fact that
every free label occurs twice. It is unchanged. The surviving quotient
is a doubled tree with total weight one after summing its leaves; here
leaf sums are precisely sum_j B_ij^2=1. Thus the same Gaussian Wick
moments and the same automorphism normalization survive, uniformly in
the chosen root and its shore.

At two fixed output roots, paired branches give

    sum_a B_ia B_ja=(B^2)_ij=:Q_ij.

If the roots are on opposite shores this sum and the paired diagrams both
vanish. At roots on the same shore, every other doubled branch is summed
away by the row-squared-norm identity. Therefore all old finite local
laws, including the actual inverse isometry U, have exactly their previous
normalization.

For own-spin marked inputs V_R, the operator-Gram proof in
`fresh_limit_injective_input_gram_2026_09_05.md`, Sections 2--6, also
reconstructs directly. At equal roots a whole isomorphic pair has weight
1+O(1/m); nonempty parity gives O(m^(-1/2)). At different roots each
pattern has an extra prefactor 1/m. If only the root-to-root parity edge
survives, consistency leaves the corresponding off-diagonal block C/m
or C^T/m, with operator norm O_L(m^(-1/2)). If another edge survives,
double the Frobenius graph. At least one single C edge remains between
free labels, and beta(C)/m^2 bounds its normalized sum. Thus the marked
input Gram is I+o_op. The output-root collision still has diagonal
covariance O(1/m), off-diagonal covariance O(1/m^2), hence covariance
operator O(1/m). Consequently

    Cov(X_T,X_U)=1_(T=U) Q+o_op.                             (4)

## 4. Global/proper cuts and nonlinear covariance

The global tree-cut proof uses only an operator bound on B and an entry
bound O(N^(-1/2)); it does not use nonzero entries everywhere. For a cut
with I crossing-isolated vertices and e internal edges, its bound becomes

    L^(crossing edges) (2m)^(I/2) m^(-e/2)=O_(T,L)(1),

because I<=e. A fixed root gains a further m^(-1/2) in every proper cut.
Collision tensor Hilbert norms are respectively O(1) globally and
O(m^(-1/2)) with the root fixed. All transported-forest proper-cut bounds
retain these exponents; N=2m changes only fixed constants.

In the nonlinear covariance proof, whole-branch pairings give the exact
old covariance factors (4). Genuinely partial matches either have two
disjoint proper-cut gains, or have one proper-cut gain multiplied by a
whole old covariance factor of Frobenius norm O(sqrt(N)). Those are the
same alternatives as in Sections 4--6 of
`resumed_bound_audit_full_nonlinear_covariance_trace_2026_09_06.md`.
Every genuinely partial contribution is O(1) in Frobenius norm and o(N)
in nuclear norm. Schur multiplication by a bounded Gram matrix is bounded
on trace class, and Q is still a correlation matrix. Thus for an old odd
response with nonlinear part R and Hermite weights w_k,

    ||Cov(R)-sum_(odd k>=3) w_k Q^(circ k)||_*=o(N).           (5)

Raw local Wick errors and marked collision errors pass through covariance
Cauchy--Schwarz exactly as before. Diagonal-free Gaussian and Rademacher
kernel covariances agree exactly. None of these steps requires same-shore
edges.

## 5. The potentially delicate old/new and full-return contractions

It is not enough to have local Gaussian moments and (5); the actual full
return has to be independent of the retained old first-chaos transports
in the ordered local limit. The two delicate diagram checks are as follows.

For Cov(V_R,P), where the source P has odd local branch count k>=3,
the only odd parity vertices are its roots a,j, and the parity degree at
j is exactly k. A free-free parity edge again gains beta(C)/m^2. If no
such edge exists, the parity graph has p length-two paths a-u-j and
possibly the edge a-j, where p=k-e>=2 and e is 0 or 1. The e=1 case is
identically zero in the bipartite array, since it contains a triangle.
In the e=0 case, both roots have the same shore and the exact sum is

    m^(-1/2) [Q^(circ p)-I],

up to a fixed normalization factor and Frobenius-small collision error.
There is no extra kernel from the doubled edges: their only residual
constraint was the consistent, already fixed coloring. The correlation
Schur bound gives operator O(m^(-1/2)) and Frobenius O(1). This supplies
the squared cross-covariance control required for the root-not-hit full
contractions, not just entrywise local convergence.

For the root-hit class, retain the exact form diag(Q J B). Its doubled
diagram has nominal O(N) size. The two P copies have separate contracted
label sets, and at most their common a-j edge can be doubled. Since k>=3,
surviving parity incidences remain at j. Any consistent diagram therefore
has a single C edge whose two labels are freely summed; inconsistent
colorings contribute zero. The beta(C)/m^2 gain makes its squared
Frobenius contribution O(sqrt(N)); additional collisions cost O(1).
Hence E||J||_F^2=o(N), and bounded left/right transports imply the required
averaged diagonal error is o(1).

These are exactly the open-contraction obligations in Sections 10--11
of the full nonlinear covariance source. Proper contractions, sign-input
Stein replacement, and marked two-spin differences use the cut bounds
above, bounded row norms, and maxentry O(N^(-1/2)); zeros are harmless.
The own-spin deterministic coefficient remains Q_ii=1. Thus the literal
old inverse K_F is unchanged, not replaced by an untested formal carrier.

## 6. Variance averaging, approximation, and the numerical value

Q=B^2 is a PSD correlation matrix even though it has two diagonal blocks.
The odd Schur theorem therefore applies without alteration:

    sum_i sqrt((B [sum w_k Q^(circ k)] B)_ii)
      >= N sqrt(sum w_k).

It permits arbitrary coherent Q and does not need a mixing Markov kernel.
Nuclear covariance errors in (5) are negligible after fixed bounded B.
The fixed-noise/fixed-softsign argument then gives the full marked return

    liminf (1/N) sum_i E H_i |(BF)_i|
      >= E H E_Z |K_F+tau Z|.                               (6)

All approximations are finite before m tends to infinity: fixed ancestor
frame and polynomial degree, fixed L, then the matrix limit, then the
Gaussian approximation and smoothing limits, and finally the countable
Gaussian closure. The exact same 21-anchor/two-Gaussian frozen policy and
rational interval certificate give the right side of (6) at least c_*.
There is no bipartite parameter in that analytic certificate.

## 7. Spectral deletion with the shores balanced

Only this special regular support is needed to remove the operator bound.
Let A_0=[[0,C],[C^T,0]]. Its ordinary bilinear norm is exactly 2 beta(C),
since the two rectangular terms can be optimized independently. The
simultaneous diagonal majorant has

    D>=A_0,-A_0,       Tr D<=2 K_G beta(C).

On a sequence beta(C)<=C_0 m^(3/2), remove coordinates with
D_ii>2 K_G C_0 sqrt(m)/epsilon. Their total number is at most epsilon*m.
If the two shores have different retained sizes, additionally remove
arbitrary vertices from the larger shore. The two retained shores then
have the common size

    k>= (1-epsilon)m.

The retained rectangular matrix C' is again a square full sign matrix.
After division by sqrt(k), its bipartite dilation has exact unit row
norm, and operator norm at most

    2 K_G C_0 / (epsilon sqrt(1-epsilon)).

Rectangular bilinear norm is monotone under row/column restriction by
zero extension in the continuous cube and independent rounding. Hence
beta(C)>=beta(C'). Apply (3)--(6) at fixed epsilon, then send epsilon to
zero. This proves (1). Boundedness along a minimizing sequence follows
already from elementary random signs; no conjecture about minimizers or
any existing Gale--Berlekamp limit is used.

## 8. Scope and checks

The proof gives a new consequence of the audited marked mechanism for
complete square bipartite sign arrays. It does not assert an extension to
all bounded-amplitude row-regular arrays, unequal rectangular aspect
ratios, sparse supports with unbounded normalized amplitudes, or a
favorable amplitude-flattening construction.

The exact finite check
`computations/decisive_audit_bipartite_tree_gram_2026_09_07.py`
checks the first nontrivial injective odd tree and its marked input,
including opposite-shore zeros, squared-factor restrictions, finite
collision counts, and the factor in (2). It is a normalization check,
not a substitute for the fixed-diagram proof above.

The reproducible run passed 64 bipartite matrices, 4384 exact marked/output
Gram entries, and 64 exact cap normalizations. The separate reconstruction
`decisive_bridge_bipartite_marked_lower_audit_2026_09_07.md` independently
checks the complete-bipartite adaptation, including its first nonlinear
tree formulas and exceptional full-return contractions.
