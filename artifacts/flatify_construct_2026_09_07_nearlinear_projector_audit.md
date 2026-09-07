# Independent audit: near-linear-rank actual-sign projector realization

2026-09-07. Adversary construction, independently reconstructed by construct. **PASS**. The canonical construction is maintained by the adversary; this note records independent checks, not a replacement proof file.

Let H be an n by n full Hadamard matrix, HH^T=nI. Let P,Q be arbitrary prescribed orthogonal projectors, with total rank r. Independently gauge H'=D_s H D_t with Rademacher row and column signs, and put

    T=(I-P)H'(I-Q).

Then T kills the prescribed subspaces exactly and ||T||op<=sqrt(n). For r log(n)/n=o(1), there exists a full-sign C with

    beta(C-T) <= O(n^(3/2) (r log(n)/n)^(1/6)+n).

All constants are universal. In particular every rank bound r<=n^(1-delta), delta>0 fixed, gives a power-saving error after absorbing the logarithm.

## Random gauge entry estimates

Keep rows i with P_ii<=mu and columns j with Q_jj<=mu. At most s<=r/mu row/column indices are discarded in total. For p>=2, Khintchine gives

    ||(PH')_ij||Lp <= C sqrt(p P_ii),
    ||(H'Q)_ij||Lp <= C sqrt(p Q_jj).

The mixed term is a decoupled degree-two Rademacher chaos. More explicitly, first condition on t and apply Khintchine in s; then apply Minkowski to the squared coefficient sum and Khintchine in t. Since every |H_kl|=1, this gives

    ||sum_kl P_ik s_k H_kl t_l Q_lj||Lp
       <= C p sqrt(P_ii Q_jj).

Thus p proportional to log(n), followed by Markov and a union bound, yields a gauge for which every kept entry satisfies

    |T_ij-H'_ij| <= epsilon
       =O(sqrt(mu log(n))+mu log(n)).

In particular no assumption about a random conditional variance is hidden in the mixed-term estimate.

## Masking, feasibility, and variance

Let T0 be T with the discarded rows and columns set to zero. Then

    beta(T-T0)<=2n sqrt(s).

Indeed split the difference into the removed-row piece and the remaining-row/removed-column piece, and use the operator norm sqrt(n), together with Boolean norms sqrt(n) and sqrt(s).

Successive orthogonal projections lose at most nr in squared Frobenius norm from H'. Every row and every column of T has squared norm at most n, so masking loses at most ns more. Consequently

    ||T0||F^2 >= n^2-n(r+s).

The entrywise mean matrix M=T0/(1+epsilon) is feasible in [-1,1]. Independent sign rounding has total variance

    V=sum_ij(1-M_ij^2)<=2epsilon n^2+n(r+s).

A Bernstein bound and union over both Boolean spin families give

    beta(C-M)<=O(sqrt(nV)+n).

Scaling costs at most epsilon n^(3/2), since beta(T0)<=n^(3/2). Combining with masking yields normalized error

    O(sqrt(epsilon)+sqrt((r+s)/n)+epsilon).

Set a=r log(n)/n and q=mu log(n). For a sufficiently small choose q=a^(2/3), so epsilon=O(q^(1/2)) and s/n<=a/q=a^(1/3). The error is O(a^(1/6)), as claimed.

## Direct-projection simplification is valid

Even without random gauges, the earlier nuclear-budget operation applies directly to T=(I-P)H(I-Q):

    Delta=H-T=PH+(I-P)HQ,
    ||Delta||* <= r sqrt(n),   ||T||op<=sqrt(n).

The generic rectangular target-contraction and rounding proof did not require T to be orthogonal. The orthogonality assumption only improved a variance constant. Therefore neither rotations nor sparse coordinate images are necessary for the earlier rank-o(sqrt(n)) operation. The random-gauge/entrywise-feasibility argument above is a substantive further improvement, not just that simplification.

## Scope

This realizes a projected Hadamard bridge as actual signs with a uniform bilinear error. It does not say that a child-energy sum plus this bridge is bounded by the favorable-flatification target. The independently proved near-orthogonal width obstruction applies to every sublinear prescribed subspace, so a separately maximized residual-operator envelope still cannot provide that payment.
