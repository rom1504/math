# Independent audit: direct projector target suffices

2026-09-07. PASS. This simplifies the arbitrary-subspace bridge construction; it does not pay its internal child energies.

Let H be an order-n full sign Hadamard, and P,Q orthogonal projections of ranks r_L,r_R. Put r=r_L+r_R and

    T=(I-P)H(I-Q), Delta=H-T.

Then ||T||op<=sqrt(n), TQ=0, PT=0, and

    rank Delta<=r,
    L=||Delta||*<=r sqrt(n).

The nuclear estimate follows by decomposing Delta=PH+(I-P)HQ. For Boolean x,y,

    |x^T T y|<=sqrt(n)||(I-P)x||||(I-Q)y||.

No localized-field error is necessary: both designated subspaces are killed exactly in the REAL target.

The already audited rectangular contraction rounding remains valid for this nonorthogonal target. Define ell_i=(sqrt(Delta Delta^T))_ii and k_j=(sqrt(Delta^T Delta))_jj, and D_ii=(1+ell_i)^(-1/2), E_jj=(1+k_j)^(-1/2). Then B=DTE is entrywise in [-1,1] and beta(T-B)<=2n sqrt(L), exactly as in the orthogonal proof.

The variance of independent mean-preserving sign recovery is

    V=n²-||DTE||F².

The nonorthogonal Frobenius deficit is paid, not omitted:

    n²-||T||F²
      =<H+T,H-T><=2sqrt(n)L.

The contraction deficit is at most 2nL because every T row and column has squared norm at most n. Hence

    V<=2sqrt(n)L+2nL<=4nL.

Union Bernstein over all 2^(2n) bilinear sign pairs therefore gives a full sign bridge C with

    beta(C-T)<=O(n sqrt(L)+n)
       <=O(n^(5/4)sqrt(r)+n).

This proves the desired actual residual bound for r=o(sqrt(n)), without any orthogonal rotations. One can alternatively compute the first deficit exactly as n(r_L+r_R)-||PHQ||F², also consistent with the estimate above.

Zero-embedding into a slightly larger Hadamard order and restricting the recovered signing works exactly as in the previously audited all-order construction. The target and rounding estimates hold on the real cube as well as its Boolean vertices by separate linearity, so zero-padding introduces no artificial Boolean padding error.
