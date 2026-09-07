# Actual low-cap signs have a bounded normalized third spectral moment

2026-09-07. Root proposed a leverage-deletion spectral uniform-integrability argument. Construct independently audits it as valid and strengthens/simplifies it by a direct Grothendieck factorization. No spectral-rounding or cap-cost conclusion is asserted.

Let A be a real symmetric hollow n by n full-sign matrix. Write

    Q(A)=max_{x in {+1,-1}^n} |x^T A x|/2,
    beta(A)=max_{x,y in {+1,-1}^n} |x^T A y|.

Let K_G denote the real bilinear Grothendieck constant. Then

    Tr |A|^3 <= K_G (n-1) beta(A)
               <= 4 K_G (n-1) Q(A).                         (1)

In particular Q(A)<=C n^(3/2) implies

    (1/n) sum_i (|lambda_i(A)|/sqrt(n))^3 <= 4 K_G C.

Thus the normalized spectral second moments are uniformly integrable, quantitatively:

    n^(-2) sum_{|lambda_i|>K sqrt(n)} lambda_i^2
       <= 4 K_G C/K.                                      (2)

This is stronger than the initially proposed O(K^(-2/3)) tail.

## Direct proof of (1)

Set |A|=(A^2)^(1/2). Every row of A has Euclidean norm sqrt(n-1), and so does every row of |A|, since

    (|A|^2)_{ii}=(A^2)_{ii}=n-1.

For n>=2, define unit vectors u_i=A_{i,*}/sqrt(n-1) and v_j=|A|_{j,*}/sqrt(n-1). Since |A| is symmetric,

    <u_i,v_j>=(A|A|)_{ij}/(n-1).

Apply the real bilinear Grothendieck inequality to coefficient matrix A and these two independently indexed families of vectors:

    sum_ij A_ij <u_i,v_j> <= K_G beta(A).

The left side is Tr(A^2|A|)/(n-1)=Tr|A|^3/(n-1). This proves the first inequality. No symmetric Grothendieck variant is used, and the nonzero diagonal of A|A| is harmless because A itself is hollow.

For completeness, beta(A)<=4Q(A): with Boolean x,y, set z=(x+y)/2 and w=(x-y)/2. Both lie in the cube, and symmetry gives x^T A y=z^T A z-w^T A w. Hollow quadratic forms attain their extrema over the cube at Boolean vertices by independent coordinate rounding. Hence each of |z^T A z| and |w^T A w| is at most 2Q(A).

The case n=1 is immediate.

## Primary theorem and normalization checked

Braverman, Makarychev, Makarychev and Naor, *The Grothendieck constant is strictly smaller than Krivine's bound*, Forum of Mathematics, Pi 1 (2013), equation (1.1), gives precisely the real matrix/unit-vector/sign formulation used above. The original research paper was read directly at

https://www.cambridge.org/core/services/aop-cambridge-core/content/view/875137460BB05E1057B9282D539AC0B6/S2050508613000048a.pdf/the-grothendieck-constant-is-strictly-smaller-than-krivines-bound.pdf

Its theorem bounds the vector objective by K_G times the bilinear sign maximum. An absolute value on that maximum is equivalent, since one whole sign family can be reversed. No current numerical estimate of K_G is needed here.

## Rank-r corollary and audit of the proposed deletion proof

For any spectral subset of size at most r, Holder gives

    sum_selected lambda_i^2
       <= r^(1/3) (Tr|A|^3)^(2/3).

Consequently its normalized energy theta satisfies

    theta <= (4 K_G C)^(2/3) (r/n)^(1/3).                 (3)

Root's original argument is also valid. Let B=P A=A P be the selected spectral part and B=L V^T with V orthonormal columns and L=A V. Each L row has norm at most sqrt(n-1), while the V row squared norms sum to r. Delete columns with leverage greater than r/(epsilon n). At most epsilon n columns are deleted. In each column,

    <A_{*,j},B_{*,j}>=||B_{*,j}||^2<=n-1,

because B is an orthogonal left projection of A. Hence the surviving inner product is at least (theta-epsilon)n^2, and its factorization has vector-norm product at most sqrt(r/epsilon). Grothendieck gives theta-epsilon<=4 K_G C sqrt(r/(epsilon n)). Taking epsilon proportional to theta proves (3), with a less sharp constant. For a threshold tail one can already improve that argument using r/n<=theta/K^2, yielding O(K^(-1)). The third-moment proof eliminates both the deletion and this self-consistency step.

## Scope and non-consequences

The theorem holds for every low-cap actual signing, without optimality, convergence, or a selectable-minimizer assumption. It gives a normalized Schatten-three bound, not just a bound on the largest eigenvalue.

It does not show that removing the high spectral tail costs o(n^(3/2)) in quadratic cap. Small normalized Frobenius energy does not supply that cap estimate. Nor does it create a full-sign replacement of a truncated matrix, preserve its diagonal or entry magnitudes, or establish favorable flatification. Covariance constructions requiring a fixed spectral bound still need to pay for truncation in the actual joint objective.

An archive text search for Schatten-three/third-spectral-moment and the specific factorization found no existing copy of this argument; this is not a claim of mathematical novelty in the literature.
