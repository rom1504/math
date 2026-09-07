# Independent audit: quenched Gaussian-sign universality

2026-09-07. Full canonical draft read:
`flatify_independent_2026_09_07_gaussian_sign_quenched_universality.md`.
**PASS**, including arbitrary deterministic offsets and the absolute parent cap.

For fixed epsilon I<=R<=KI with unit diagonal, actual Gaussian signs C=sign(G), and Gaussian Y with exactly Cov(Y)=Cov(C), the canonical theorem proves

    |E Phi(C)-E Phi(Y)|
       <=C_(epsilon,K)|lambda|^3 d log(d+1)^(3/2),

uniformly over the number of configurations, offsets, and coordinate features bounded by one. This is a quenched pressure comparison, not the earlier fixed-pair MGF or its union bound.

## Independent checks of the substantive steps

1. Conditional on W in G=W+sqrt(epsilon)Z, the sign coordinates are independent. Replacing them by f(W_i)+sqrt(1-f(W_i)^2)xi_i matches their first two conditional moments. The pressure third partial derivatives are uniformly O(|lambda|^3), so conditional Lindeberg costs O(|lambda|^3 d), irrespective of the configuration count. The resulting F has the exact original sign covariance, including its diagonal.

2. The local smooth function g(w,xi)=f(w)+sigma(w)xi has both gradient and Hessian bounded by C_epsilon(1+|xi|). The scalar sigma derivatives really are bounded: inverse Mills inequalities bound its logarithmic derivatives polynomially, while sigma itself has Gaussian decay. The proof's r'=r(r-z)<=2 for z>=1 follows from z<=r<=z+1/z. No unproved uniform smoothness at the tails is being assumed.

3. The block covariance T=diag(R-epsilon I,I) in the grouped coordinates has uniformly bounded operator norm. OU smoothing of a function of block i stays a function of block i, even when the blocks are correlated. Thus the Stein kernel Gamma_ij=a_i^T T_ij b_j has derivative support only on blocks i,j. The integration-by-parts orientation in the canonical equation (6) is correct. Neither symmetry of Gamma nor invertibility of T is required.

4. In the covariance correction to the pressure interpolation, differentiating Gamma produces exactly two star centers: coefficients T_ij T_ik or T_ij T_jk. The derivative of the pressure Hessian is the third Gibbs cumulant. Its expansion has five expectation products with total absolute coefficient sum six, representable by at most three replicas. Once those replicas and all environments are fixed, each term is a product of one i-, one j-, and one k-indexed scalar. This remains true for repeated replica labels, coinciding indices, and different OU environments.

5. Consequently the essential contraction is pointwise:

       |sum_i a_i (B b)_i (D c)_i|
          <=||a||infinity ||B||op ||D||op ||b||2 ||c||2.

Each block-coordinate submatrix B,D of T has norm at most ||T||op. With bounded feature coordinates, this is at most d||T||op^2 times the product of three local derivative maxima. Dependence of the replica spins on every Gaussian coordinate does not affect this deterministic inequality. Thus the proof does not replace an operator norm by an absolute-entry norm or silently assume spin/coefficient independence.

6. No auxiliary Gaussian truncation is needed. The derivative product is bounded by C_epsilon(1+max|xi|+max|xi_s|)^3, whose expectation is O(log(d+1)^(3/2)) uniformly in OU time. Both arrays have standard Gaussian marginals; their dependence is harmless. Pressure derivatives are bounded, and all remaining functions have polynomial Gaussian growth, justifying integration by parts and dominated limits.

## Ground-state rate

With d=n^2, lambda=beta/sqrt(n), and absolute-parent configurations (sigma,x,y), the normalized expected cap error is

    O(1/beta)+C_(epsilon,K) beta^2 n^(-1/2) log(n+1)^(3/2).

Taking

    beta=n^(1/6)/sqrt(log(n+1))

balances the two terms and proves the quantitative rate

    O_(epsilon,K)(n^(-1/6) sqrt(log(n+1))).

The canonical draft's simpler beta=n^(1/12) also proves convergence; the displayed choice sharpens it. No upper bound on beta is required by the proof, and this choice still has lambda tending to zero.

## Exact remaining scope

The theorem transfers the expected absolute parent maximum, with the unchanged child energies present as arbitrary deterministic offsets. It therefore genuinely escapes the central-shell obstruction to separate fixed-pair Chernoff/union bounds.

It does not evaluate or upper-bound the covariance-matched Gaussian parent maximum at the favorable-flatification target. Constants require fixed epsilon,K. For opposite-spectral covariance R=I-rho A tensor D/(||A||op||D||op), any fixed rho<1 meets those constants, but a large child operator denominator weakens the energy-dependent covariance. The new spectral fourth-moment theorem does not by itself remove that loss or settle the Gaussian optimization.
