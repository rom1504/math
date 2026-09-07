# Independent audit of growing-rank residual energy width

2026-09-07. Root's
`flatify_director_subspace_residual_width_2026_09_07.md` passes independent
reconstruction in full, including its uniform rank range and power errors.

The elementary operator bound is particularly useful: each row a_i of A
lies in the cube, so ||(A^2)_i||1=max_s a_i^T A s<=beta(A). Since A^2 is
symmetric, its operator norm is bounded by this maximum row sum. Thus
||A||op^2<=beta(A)<=4Q(A), without interpolation or an eigenvector loss.

For the split k=floor(n/3), l=n-k, the Gaussian signs sign(g_i) and
sign((Bg)_j) have correlation B_ji/sqrt(k). Multiplying by the actual
edge sign gives the same positive contribution for every cross edge.
Hence E L=kl*(2/pi)*arcsin(1/sqrt(k)), with leading coefficient
4/(3pi sqrt(3)). The Gaussian pre-sign covariance has unit diagonal and
operator norm 1+||B||op^2/k=O_C(sqrt(n)).

The entrywise arcsine covariance estimate is rigorous: its positive
odd-power coefficients sum to one after normalization. Schur multiplication
by a correlation matrix is positive and unital and therefore contracts
the operator norm on symmetric matrices. Repeated Schur products have
norm at most the original correlation norm. Perfectly correlated rows
and endpoint entries +/-1 cause no issue because the positive coefficient
series converges at one.

Thus BOTH witnesses have expected squared projection O_C(r sqrt(n)).
The union Markov bound for threshold tau*n is O_C(r/(tau sqrt(n))).
The cross energy need not be nonnegative pointwise; nevertheless it obeys
|L|<=Q(A), because the two full energies differ by 2L. Therefore discarding
bad projections loses at most Q(A) times their probability from E L.
This is sufficient for an actual good sample with the claimed half-width.

For r=o(sqrt(n)), choose tau->0 more slowly than r/sqrt(n). For
r<=n^(1/2-2delta), tau=n^-delta gives the stated O_C(n^-delta+n^-1)
relative error. The covariance and expectation estimates are independent
of the chosen subspace and uniform in the cap constant.

Finally, half-width rather than one-sided energy is correctly used when
combining arbitrary child polarities. The low-projection sets retain
signed-energy interval widths, while the scalar residual upper envelope
is at least 1-tau on every pair in those sets. This forces the stated
certificate lower value 1+8/(3pi sqrt(3))-o(1). It is NOT a lower bound
on the actual bridge energy at those witnesses. Joint bridge-child
alignment, and operations changing internal edges, remain outside this
no-go scope.
