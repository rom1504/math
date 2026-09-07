# Bounded-entry matrices: Schatten-four versus bilinear sign norm

2026-09-07. Root asked whether actual low cap forces a bounded normalized fourth spectral moment. **Yes.** The following general inequality follows from Grothendieck factorization and the real bilinear Grothendieck inequality:

    ||A||_{S4}^2 <= K_G ||A||_max beta(A),                 (1)

for every real rectangular matrix A. Here beta(A)=max_{x,y signs}|x^T A y|, ||A||_max=max_ij |A_ij|, and K_G is the real bilinear Grothendieck constant.

In particular, for a hollow symmetric full signing with Q(A)<=C n^(3/2),

    Tr A^4 <= K_G^2 beta(A)^2 <= 16 K_G^2 C^2 n^3.        (2)

This strengthens the third-moment theorem in `flatify_construct_2026_09_07_spectral_third_moment.md`.

## Factorization proof

Let G(A) be the vector relaxation of beta(A): the supremum of sum A_ij <u_i,v_j> over unit vectors. The real bilinear Grothendieck inequality is G(A)<=K_G beta(A).

The corresponding factorization gives probability vectors p,q and a matrix T such that

    A=D_p^(1/2) T D_q^(1/2),     ||T||op<=G(A).           (3)

Zero weights cause no problem: the corresponding rows or columns of A vanish and can be omitted. For completeness, (3) follows directly from finite SDP duality. The dual of the vector relaxation minimizes (sum u_i+sum v_j)/2 subject to the block matrix

    [ D_u   -A  ]
    [ -A^T  D_v ] >= 0.

Rescale u by a positive scalar and v by its reciprocal to equalize their sums; this preserves feasibility by block-diagonal congruence and can only decrease the objective. At the optimum both sums are G(A). The Schur-complement/contraction characterization of the block positivity then gives (3). An approximation with strictly positive weights gives the same result if needed. Strong SDP duality holds because the dual has a strictly feasible point, obtained with sufficiently large scalar diagonal blocks.

Write M=||A||_max. Using (3) for the middle factor,

    AA^T A=(A D_q^(1/2)) T^T (D_p^(1/2) A).

Every row of A D_q^(1/2) has norm at most M because sum q_j=1. Every column of D_p^(1/2) A has norm at most M because sum p_i=1. Absorb T^T into one vector family. Thus AA^T A has a vector factorization with product of the two maximum vector norms at most G(A) M^2.

Apply the definition of G(A), now to coefficient matrix A and this factorization:

    <A,AA^T A> <= G(A)^2 M^2.

The left side is Tr[(A^T A)^2]=||A||_{S4}^4. Taking square roots and using G(A)<=K_G beta(A) proves (1). Equivalently this is two applications of Grothendieck, but retaining G(A) makes the constants transparent.

## Normalization and source

The real bilinear matrix/unit-vector/sign theorem used is equation (1.1) of Braverman, Makarychev, Makarychev and Naor, *The Grothendieck constant is strictly smaller than Krivine's bound*, Forum of Mathematics, Pi 1 (2013), read directly at

https://www.cambridge.org/core/services/aop-cambridge-core/content/view/875137460BB05E1057B9282D539AC0B6/S2050508613000048a.pdf/the-grothendieck-constant-is-strictly-smaller-than-krivines-bound.pdf

The factorization above was derived explicitly from finite SDP duality, so no unverified variant of a factorization theorem is required. The equivalent probability-measure factorization is also discussed by Pisier, *Grothendieck's theorem, past and present*, https://webusers.imj-prg.fr/~gilles.pisier/grothendieck.pdf . No numerical constant estimate is needed.

For hollow symmetric A, beta(A)<=4Q(A) by polarization on the cube, as proved in the third-moment artifact. This proves (2).

## Consequences and strict scope

Under Q(A)<=C n^(3/2), let D=16 K_G^2 C^2. Then

    n^(-2) sum_{|lambda|>K sqrt(n)} lambda^2 <= D/K^2,
    #{|lambda|>K sqrt(n)} <= D n/K^4,
    n^(-2) sum_top_r lambda^2 <= sqrt(D r/n).

These statements hold for all actual low-cap matrices, not merely minimizers. They imply spectral uniform integrability with a quadratic tail rate. They do not imply that spectral truncation has small quadratic cap, is entrywise feasible, or can be rounded at a favorable cap cost. No convergence or favorable-flatification conclusion is claimed.

The inequality is scale-consistent for arbitrary matrices: scaling A scales both sides of (1) quadratically because the maximum-entry factor is included. An archive search did not locate this exact argument; no claim of literature novelty is made.

## Fourth moment is the endpoint for general near-minimizers

Start from any exact minimizing full signing A_n. Choose s=floor(n^(3/4)/log n) vertices and change all internal edges of this set to +1, obtaining C_n. At most binom(s,2) edges change and each changed coefficient changes the quadratic form by at most 2. Therefore

    M_n <= Q(C_n) <= M_n+s(s-1)=M_n+o(n^(3/2)).

The unit vector constant on the planted set and zero elsewhere has Rayleigh quotient s-1, so lambda_max(C_n)>=s-1. For each fixed p>4,

    n^(-1-p/2) Tr|C_n|^p
       >= n^(-1-p/2)(s-1)^p
       ~ n^(p/4-1)/(log n)^p -> infinity.

Hence no uniform normalized p-th spectral moment bound for p>4 follows from additive near-minimality. This does not establish such a counterexample among exact minimizers; that narrower class is not settled by planting. It also shows why the O(n^(3/4)) operator bound cannot be improved by a fixed power for arbitrary near-minimizers.
