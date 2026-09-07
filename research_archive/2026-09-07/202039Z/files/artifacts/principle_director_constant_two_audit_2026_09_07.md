# Independent reconstruction: the constant-two independent-frame floor

2026-09-07. **PASS after reconstruction, not by inherited audit verdict.**
Root read the full canonical constant-two proof and independently derived
the following chain.

1. For C_ij=S_ij eta_ji and D_ij=S_ij eta_ij eta_ji, conditioning on D
   leaves one fair orientation per unordered edge and gives
   C_ji=S_ij D_ij C_ij. Choosing a positive maximizing y of D and
   x_i=sign((Cy)_i), with fair ties, gives exactly
   E H_S(x)=lambda_n^2 H_D(y), not a loss of Q(S).
   The shared-edge field calculation uses disjoint remaining random edges.
2. The Gaussian covariance p(D/sqrt(n))^2+epsilon I is positive definite.
   Root reconstructed the four-open-walk count: a tree union traverses the
   i--j path at least four times, while a cyclic union has vertices<=edges.
   Both give bounded total off-diagonal fourth moment. Diagonal variance,
   semicircle trace and Frobenius normalization errors follow at fixed degree.
   The arcsine remainder is summable by second/fourth moment Cauchy--Schwarz.
   The beta integral for p=(x+2)^r gives2-6/(2r+3). Thus the one-sided
   iid-child value is at least2/pi in the ordered limits n,epsilon,r.
3. The old child plus opposite-child witness and optimized bridge therefore
   yield (2/pi+sqrt(2/pi))/(2sqrt2). Each directed port changes two actual
   parent edges, giving the claimed bounded-differences high-probability law.
4. Uniform majority repair has the exact uniform balanced-row marginal.
   Its conditional mean is (1-p_i)eta_i-a_i and its covariance is bounded
   by8p_i I. Root multiplied these means directly to recover both C and D
   formulas, including the terms containing the unchanged seed S.
5. For D, exposure of repaired row i produces an actual symmetric star:
   earlier neighboring rows are realized, later ones are conditional means.
   The predictable covariance bound therefore includes all cross terms.
   Conditional on the original good event, increments and total quadratic
   variation are deterministic bounds. Root read Tropp's primary Theorem1.2
   and rectangular Corollary1.3 at https://arxiv.org/pdf/1101.3039 and checked
   dimensions, both tails, and variance conventions. They give
   O(n^(1/4)log^(5/4)n) operator repair, hence o(n^(3/2)) cap change.

The theorem applies to every child selected BEFORE independent frame
randomness, including exact minimizers; it is not simultaneous over seeds
chosen from the frames. It rules out the TYPICAL independent two-frame law,
not rare correlated completions. There is no assertion of a half-floor
for the already constructed correlated anti-invariant subhalf family.

The quantitative result is a real parent cap obstruction above .50717,
not an annealed-pressure interpretation. It leaves the original convergence
problem and arbitrary global rewriting open.
