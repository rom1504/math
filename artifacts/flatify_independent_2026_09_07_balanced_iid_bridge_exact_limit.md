# Exact iid balanced-bridge limit from a checked primary theorem

2026-09-07. Imported ensemble theorem, not convergence of the original minima.

Primary source: Hong-Bin Chen, Victor Issa, Jean-Christophe Mourrat,
[Free energy of non-convex multi-species spin glasses with centered Ising
spins](https://arxiv.org/html/2606.16636v1), June 2026 preprint, Theorem 1.1,
Definition 9.1 and Proposition 9.3. The generic theorem was already catalogued
in `retrieval_panel_2026_08/spin_glass_toolkit.md`, section 15. The explicit
normalization below is its application to the current iid-bridge mechanism.

Let B_m have independent standard Gaussian entries, and set

    X_m(x,y)=m^(-1/2) sum_{i,j<=m} B_ij x_i y_j.

For N=2m, the paper's species overlaps are
a_1=N^(-1)sum_i x_i x'_i and a_2=N^(-1)sum_j y_j y'_j.
Then Cov(X_m(x,y),X_m(x',y'))=N xi(a), where xi(a)=2a_1 a_2.
The proportions are lambda=(1/2,1/2). Its balanced comparison is

    xi(r lambda)=r^2/2=:xi_*(r),
    2a_1 a_2 <= a_1^2+a_2^2
              =sum_s lambda_s xi_*(a_s/lambda_s).

All hypotheses are satisfied: centered symmetric Ising single-site measures,
positive fixed proportions, a centered Gaussian Hamiltonian, polynomial
covariance, and no external field. Proposition 9.3 therefore identifies the
limiting pressure with that of the ordinary SK model having covariance
N r^2/2. The paper uses a negative, probability-normalized, variance-corrected
pressure at sqrt(2t)=beta. Both models have xi(lambda)=1/2, so these identical
normalization terms cancel in the comparison.

Write P_SK for the limiting positive SK ground-state energy per spin in this
normalization. Uniformly in N, ordinary log-sum-exp differs from the positive
maximum by at most N log(2)/beta. Taking N to infinity at fixed beta, then
beta to infinity, proves

    E max_{x,y} X_m(x,y)/(2m) -> P_SK.

The bipartite maximum equals its absolute maximum by y -> -y. Hence

    E max_{x,y}|sum B_ij x_i y_j| / (2m)^(3/2)
        -> P_SK/sqrt(2).

The same conclusion holds for iid unbiased Bernoulli signs. At each FIXED
beta, the previously checked softmax Lindeberg estimate is O(beta^2/sqrt(m))
per spin. Transfer the finite-temperature limit first, then send beta to
infinity; no unproved uniform-in-temperature convergence rate is used.
Bounded differences (or Gaussian concentration) turns convergence of the
mean into convergence in probability. Thus this is an exact iid ensemble
limit, whereas the earlier Slepian argument alone supplied only its lower
bound.

The explicit martingale certificate in
`flatify_independent_2026_09_07_heat_martingale_sk_bound.md` gives

    P_SK/sqrt(2) >= Gamma(3/4)/(sqrt(2*pi) Gamma(5/4))
                 =0.539352601188379356... .

For any within-block completions, even chosen after seeing B, global reversal
of one block changes the bridge energy's sign and leaves the internal energy
unchanged. Therefore the full original absolute cap is at least this bridge
maximum. Adaptive child selection cannot remove the obstruction.

## Finite-species row-regular extension

More generally, fix positive species proportions lambda_s and a symmetric
nonnegative variance matrix v_st satisfying sum_t v_st lambda_t=1 for every
s. Independent Gaussian edges with these variances and ambient normalization
sqrt(N) have covariance N xi(a)+O(1), where

    xi(a)=1/2 sum_{s,t} v_st a_s a_t.

The O(1) diagonal correction is a spin-independent Gaussian and does not
affect the expected positive maximum or expected log partition function.
Writing a_s=lambda_s r_s and applying 2r_s r_t<=r_s^2+r_t^2 gives

    xi(r lambda)=r^2/2,
    xi(a)<=1/2 sum_s lambda_s r_s^2.

Thus the same theorem gives the SK limiting pressure and ground energy for
every such fixed finite-species row-regular variance profile, including
profiles with zero cross variances. Independent Bernoulli signs with the
corresponding bounded amplitudes have the same limit by fixed-temperature
Lindeberg. This is a theorem about independent disorder, not about optimized
signs. In particular it does not control the signed drift of our optimized
partition-reveal pressure.

The imported result is a recent preprint; this note checks its stated
hypotheses and deductions, not its entire long proof.
