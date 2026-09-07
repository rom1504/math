# Universal actual-sign instability: a fourth-moment variance theorem

2026-09-07. **Proved; independently reconstructed by director and invent.** This complements, rather than replaces, the exponential lower-tail theorem in `principle_construct_2026_09_07_universal_stability_audit.md`.

Let A be any symmetric hollow full-sign matrix of order n. For a uniform Boolean word x define its total negative local field

    I_A(x)=sum_i [-x_i(Ax)_i]_+.

There is an absolute constant K such that

    Var I_A <= K [n^2+Tr(A^4)/n].                          (1)

Consequently, if Q(A)<=C n^(3/2),

    || I_A/n^(3/2)-1/sqrt(2pi) ||_(L2)
       <=K'(1+C)n^(-1/2).                                (2)

The constants are uniform over actual matrices. More generally the right side in (2) can be replaced by K'[n^(-1/2)+Q(A)/n^2]. Thus Q(A)=o(n^2) already implies the same universal L2 law, even without a bounded n^(3/2) cap coefficient. No Gaussian replacement of the actual matrix is involved.

## 1. Convex cube Poincare eliminates nonsmooth boundary errors

Put F(x)=||Ax||_1 and s(x)=sign(Ax), with sign(0)=0. The vector A^T s(x) is a legitimate convex subgradient of F at x. If x^(j) reverses coordinate j, convexity gives

    F(x)-F(x^(j)) <=2x_j (A^T s(x))_j.

Let Delta_j be the left side. Reversal symmetry of the uniform cube gives E Delta_j^2=2E(Delta_j)_+^2. The cube Poincare inequality, proved directly by tensorization or the Walsh expansion, therefore yields

    Var F <=(1/4)sum_j E Delta_j^2
           <=2E ||A^T s(x)||_2^2.                         (3)

This does not differentiate sign across a zero or approximate a discrete derivative by a smooth one. In particular, no separate small-ball boundary error is accumulated in (3).

Set K_s=E[s(x)s(x)^T]. It is positive semidefinite. Since A is symmetric,

    E ||A^T s(x)||_2^2 =Tr(A^2 K_s).                      (4)

## 2. Uniform elementary sign-covariance bound

For distinct i,j,

    |(K_s)_ij| <=K_0 [ |(A^2)_ij|/n+n^(-1/2) ].          (5)

Here is a direct scalar proof, extending the archived two-binomial argument in `flatify_independent_2026_09_07_rademacher_residual_width.md` to hollow rows.

Remove coordinates i,j from the two row sums. The remaining sums X,Y have the same q=n-2 sign-coefficient locations. Each full row differs from its remaining sum by one sign of magnitude one. A sign can change under this addition only when the remaining sum has absolute value at most one. The central binomial bound makes each such event O(q^(-1/2)). It is therefore enough to estimate E sign(X)sign(Y), up to that error.

Gauge the coefficients of X to +1. Let a coefficients of Y agree and b disagree, so a+b=q and a-b=(A^2)_ij. Then

    X=U+V, Y=U-V,

where U and V are independent fair-sign sums of lengths a,b. If |a-b|>=q/2, the trivial correlation bound one already proves (5). Otherwise a,b>=q/4.

The scalar iid Berry--Esseen inequality gives uniform CDF errors O(q^(-1/2)) for U,V compared with independent Gaussian variables of variances a,b. Passing to absolute values costs at most a factor two. Replacing the two marginal distributions successively in Pr(|U|>|V|) costs O(q^(-1/2)): the outer Gaussian absolute-value CDF is a bounded monotone function, so integration against a distribution with small Kolmogorov distance preserves that order of error. Boundary ties have probability at most Pr(X=0)+Pr(Y=0)=O(q^(-1/2)).

Since sign(X)sign(Y)=sign(U^2-V^2), the Gaussian answer is

    (2/pi) arcsin((a-b)/q),

by the elementary planar Gaussian angle identity. Arcsine is Lipschitz on [-1/2,1/2]. This proves (5), after replacing q by n with harmless absolute factors. Finitely many small n are absorbed into the constant.

For a primary quantitative source, I independently read [Shevtsova's scalar iid result](https://arxiv.org/pdf/1111.6554), Theorems 1--2 and Corollary 1. Fair signs have mean zero, variance one, and third absolute moment one, so its uniform CDF error is an absolute constant times the inverse square root of the number of summands. Only this scalar statement is used here, not multivariate Berry--Esseen or covariance approximations for absolute fields.

Squaring (5), summing, and treating the diagonal by |(K_s)_ii|<=1 gives

    ||K_s||_F^2 <=K_1 [Tr(A^4)/n^2+n].                    (6)

## 3. Variance and the exact universal mean

Write T=Tr(A^4). Combining (3), (4), (6), and Frobenius Cauchy--Schwarz,

    Var F <=K_2 sqrt(T) sqrt(T/n^2+n)
           <=K_3 [T/n+n^2].                             (7)

The final inequality follows from sqrt(nT)<= (n^2+T/n)/2.

The identity

    I_A(x)=[F(x)-x^T A x]/2

is exact. Hollowness and independence give Var(x^T A x)=2n(n-1). The variance inequality for a difference and (7) now prove (1).

Also, independently of A,

    E I_A = (n/2) E|S_(n-1)|
          = n^(3/2)/sqrt(2pi)+O(sqrt(n)),                (8)

where S_l is a sum of l fair signs. The first equality uses E(x^T A x)=0 and the exact common row-sum distribution. The error follows from the central-binomial formula and Stirling's estimate.

The previously proved actual-sign Schatten inequality gives

    Tr(A^4) <=K_G^2 beta(A)^2 <=16K_G^2 Q(A)^2.

Substituting it into (1) and using (8) yields both (2) and the broader bound stated above.

## 4. What is and is not supplied

This is a quantitative law of large numbers for total instability under UNIFORM spins. It is not a Gibbs-state universality theorem, a cap upper bound, or convergence of M_n/n^(3/2).

Its Chebyshev tails are polynomial. The director's separate diagonal-majorant/Talagrand argument proves exponentially small lower tails below every fixed level v<1/sqrt(2pi) when Q(A)<=Cn^(3/2); that stronger tail result retains independent content.

The proof reuses the established actual S4 theorem and the existing scalar two-binomial sign-correlation estimate. The new step is their combination with the exact convex cube Poincare estimate (3), giving the variance scale n^2 without an absolute-field multivariate CLT or a loss from nonsmooth sign boundaries.
