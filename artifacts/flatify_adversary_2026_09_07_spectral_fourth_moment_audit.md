# Independent audit: fourth spectral moment and residual-width upgrade

2026-09-07. PASS. Audited `flatify_construct_2026_09_07_spectral_fourth_moment.md` and independently read the primary real bilinear Grothendieck statement, equation (1.1), in [Braverman--Makarychev--Makarychev--Naor](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/875137460BB05E1057B9282D539AC0B6/S2050508613000048a.pdf/the-grothendieck-constant-is-strictly-smaller-than-krivines-bound.pdf).

Let G(A) be the unit-vector relaxation of beta(A). The SDP objective has off-diagonal coefficient blocks A/2 and A^T/2, so its dual is exactly half the sum of diagonal vectors u,v subject to [[D_u,-A],[-A^T,D_v]] positive semidefinite. Reciprocal rescaling of the two diagonal blocks is a congruence preserving feasibility; balancing their sums gives both sums G(A). Normalizing u,v to probability vectors therefore yields

    A=D_p^(1/2) T D_q^(1/2), ||T||op<=G(A).

This is not a factorization imported with a hidden constant. Strict feasible points give the required finite-dimensional strong duality. Zero weights force the corresponding A row or column to vanish; A=0 is trivial.

For M=max|A_ij|, every row of A D_q^(1/2) and every column of D_p^(1/2)A has Euclidean norm at most M. The exact identity

    AA^T A=(A D_q^(1/2)) T^T (D_p^(1/2)A)

thus represents AA^T A as inner products with product of maximum vector norms at most G(A)M². Pairing these vectors with the ORIGINAL coefficient matrix A costs at most another G(A). Vectors of norm at most one may be extended to unit vectors without changing cross inner products, so there is no unit-ball/unit-sphere issue. It follows that

    ||A||S4^4=<A,AA^T A><=G(A)²M²
       <=K_G² beta(A)²M².

Dimensions, transposes, and all factors of two check. This is neither circular nor dependent on optimizing-child regularity. For hollow signs with Q<=C n^(3/2), it proves trA^4<=16K_G²C² n³.

The stated spectral tail and top-r energy bounds follow directly from this fourth moment. The planted clique of size n^(3/4)/log n changes cap by o(n^(3/2)) yet violates every normalized p-th moment bound for fixed p>4; its Rayleigh quotient is exactly s-1 on the constant supported vector. This counterexample concerns near-minimizers, not exact minimizers.

## Stronger consequence for the Rademacher residual-width theorem

In `flatify_independent_2026_09_07_rademacher_residual_width.md`, equation (8) now gives

    ||E ss^T||F²<=C[q^(-2)tr A^4+|P|²/q+|P|]=O(n).

The mixed and original-seed covariance blocks also have squared Frobenius norm O(n), so the whole witness covariance satisfies ||K||F=O(sqrt(n)). For any rank-r projection, both witness projection energies have expectation at most C sqrt(rn).

Hence the same deletion argument proves half-width at least

    [(2/3)sqrt(2/(3pi))-o(1)] n^(3/2)

near the orthogonal complement of EVERY prescribed r=o(n) subspace. No new CLT estimate is needed. In particular the separately paid residual-norm envelope has the same 1.6142118... lower coefficient throughout all sublinear ranks. This is still only a scalar-certificate obstruction, not a lower bound on an actual selected bridge at those witnesses.
