# Director reconstruction of the full nonlinear response channel

Date: 2026-09-06. This reconstructs the new argument in
`resumed_bound_audit_full_nonlinear_covariance_trace_2026_09_06.md`.
It does not rely on the response researcher's audit verdict. The original
convergence problem remains open. Numerical consequences require their own
exact arithmetic certificates.

## Covariance module

1. The branch-matching graph records actual slot contractions, with all
   branch marked degrees odd. A component that is neither a star nor an
   edge contains a P4. Contracting its two disjoint outer edges first gains
   two proper fixed-root flattenings, hence O(1/n). Further tensor-network
   contraction is bounded by Hilbert norms; no dimension factor is inserted.
2. A lone nontrivial star has an odd number of leaves, at least three. Since
   both response monomials have local degree at least three, at least two
   other whole-branch edges remain. One bounded-operator old Gram factor
   has Frobenius norm O(sqrt(n)); multiplying it by the star's uniform
   O(n^(-1/2)) entry bound gives O(1) Frobenius norm. This argument would
   fail for first local chaos, which is correctly excluded.
3. Whole-branch terms use the actual old covariance-operator theorem, not
   merely joint Gaussian moments at two roots. Cross-type Gram factors
   have bounded Schur multiplier norm by their row-vector factorizations.
   The normalized Hermite factors cancel the branch pairing multiplicities.
4. The conclusion is normalized NUCLEAR convergence, not o(1) operator
   convergence. Averaged local L2 errors therefore suffice to transfer
   covariance: ||E UV^T||_*<=sqrt(E||U||^2 E||V||^2). Deleting repeated
   marked slots costs vanishing averaged Hilbert norm, after which Gaussian
   and Rademacher covariances agree exactly.
5. The scalar inequality |sqrt(a)-sqrt(b)|<=sqrt(|a-b|), summed against
   rows of B, transfers the proved Schur mean-standard-deviation lower
   bound. No operator-Lipschitz square-root claim is used.

## Independent innovation and the exceptional first chaos

6. Write each old output X_T=B V_T+E_T, using the exact own-marked even-root
   input and its small covariance-operator collision error. In Cov(V_T,P),
   the new forest root j is excluded from BOTH spin sets and has parity
   degree exactly k>=3. A free-free parity edge yields the O(1) Frobenius
   estimate directly. Without one, the graph is an optional aj edge and
   p=k-e>=2 length-two paths. The resulting normalized Q Schur powers
   or B circ Q powers have O(n^(-1/2)) operator norm. Thus all-pairs
   Cov(BP,X_T) has bounded Frobenius norm.
7. Consequently same-root covariances vanish in average, so the Gaussianized
   higher-chaos channel is independent of the local old family in the
   averaged comparison. This removes, rather than assumes away, the possible
   old linear drift left open by the first forest flattening theorem.
8. B X_T itself may still contain coherent own-spin atoms and is not declared
   Gaussian. Proper mixed contractions into it are small by the new channel's
   own proper flattenings. The full contraction must be split separately.
9. For the root-not-hit part, first replace the even child feature by its
   pure top input-chaos kernel. The positive covariance-square identity
   then bounds the full contraction by squared old/new cross-covariances
   plus an averaged vanishing error uniform in the other root. The all-pairs
   Frobenius estimate in item 6, bounded Q entries, and sum Q_ia^2=O(n)
   suffice after averaging. The own-free Gaussian creation inequality
   handles the remaining sum of Q_ia N_a without an n-fold absolute loss.
10. The root-hit part is diag(Q J B). My separate complete counting proof in
    `resumed_full_forest_root_hit_director_2026_09_06.md` gives
    E||J||_F^2=O(sqrt(n)), and hence averaged root-hit square O(n^(-1/2)).
    The essential fact is that at least four incidences at the new root
    survive parity in a leading doubled diagram. Residual marks at that
    root are lower label patterns, not missing leading terms.
11. These are the unsymmetrized contractions needed for block polynomial
    moment independence. In a full-channel contraction the surviving slots
    are symmetric already. No joint Gaussian assertion about the exceptional
    first-chaos background, or moment determinacy of that background, is used.

## Signs, bounded feasible functions, and limit order

12. Squarefree forest deletion is an orthogonal projection on the common
    marked tensor space. Thus its global root-map norm stays bounded, even
    though proper flattened norms may only remain small on average. The
    averaged maximum influences tend to zero; total influences stay bounded.
13. In unmarked hybrid replacement Y=B X_T is squarefree and appears only
    linearly. Fourth derivatives have only Y D^4f and 4(DY)D^3f. With local
    derivative scales c_(i,a), sum_a c_(i,a)^2=O(1), and average
    max_a c_(i,a)^2=o(1), both terms vanish after summation. This does not
    require small individual influences of Y.
14. In marked tests the direct two-spin identity retains the own sign exactly.
    Its off-diagonal sum is controlled by the old derivative n^(-1/2),
    the Q row l2 norm, and the channel's total influence. The diagonal
    gives the creation adjoint K_F. Root-spin removal is only applied to
    the new low-influence channel, not to the coherent first-chaos term.
15. At a fixed odd polynomial P, the FOREST MAIN representative of the
    new channel B(P-P1P) has uniformly bounded variances and the proved
    old independence. The raw channel is close only in averaged L2;
    normalized nuclear covariance does not imply uniformly bounded raw
    root variances. Use the forest main in the auxiliary-noise density
    domination, and transfer by averaged L2 at each fixed softsign scale.
    Its mean standard deviation is at least ||P-P1P||_2-o(1), since scalar
    square-root comparison transfers this from the raw variance result.
    The restricted-test/softsign proof
    then lower-bounds average E H|BP| by E H E|U*P+||P-P1P||_2 N|.
    The auxiliary Gaussian smoothing still comes BEFORE bounded-test
    approximation, including at roots of zero variance.
16. For actual bounded feasible F,H, only the exact means +/-F+H sign(BF)
    are used as cube points. An odd polynomial P need not be feasible:
    it is an analysis approximation, and fixed ||B||op bounds
    |average E H|BF|-average E H|BP|| by L||F-P||_(2,n).
    After the fixed-P matrix limit, improve P in Gaussian L2. Both U*P
    and the nonlinear norm converge by orthogonal-projection contractions.
17. The resulting Gaussian functional has no dependence on L. Principal
    deletion therefore removes the fixed operator bound last. Infinite
    old-coordinate feasible pairs can afterward be conditioned and smoothed
    jointly, preserving feasibility, just as in the all-mask escape proof.

## Conclusion of this reconstruction

The full-channel mechanism passes these checks. In particular it supports
the stronger universal functional

    liminf M_n/n^(3/2) >= E H E_N|U*F+tau_F N|,
    tau_F^2=||F||_2^2-||P1F||_2^2,

for fixed old-coordinate odd bounded F and even H with |F|+H<=1, and by
ordered finite approximation for the Gaussian closure. The new noise is
the entire nonlinear local Hermite residual, not merely an edge projection.
A separate independent researcher completed its reconstruction in
`resumed_response_full_channel_independent_audit_2026_09_06.md`. The exact
decimal has its own numerical audit. This one-step theorem does not
refresh the now-dependent center or provide an all-order upper construction.
