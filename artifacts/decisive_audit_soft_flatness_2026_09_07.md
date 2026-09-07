# Independent audit of the soft-flatness variational bridge

Date: 2026-09-07. Full mathematical reconstruction PASS for
`decisive_independent_soft_flatness_variational_bridge_2026_09_07.md`.
No original convergence claim is implied without its stated new limit
obligation.

The strong-variance bound is legitimate for arbitrary independent centered
finite-variance coordinates, even without bounded tails. For the smooth
maximum F_lambda, the global bound partial_e^2 F_lambda<=lambda permits
Taylor replacement of coordinate e by zero with expected cost at most
lambda Var(Z_e)/2. The first-order term is independent of Z_e and has
zero expectation. Summing these replacements and optimizing lambda gives
sqrt(2 log|V| sum Var Z_e). Thus no Bernstein remainder is needed.

Using at most 2^n projectively reduced signed quadratic states and
sum Var(A_e-B_e)=d delta(B) gives exactly

    m_n<=q_n(B)+sqrt(a_n delta(B)),
    a_n=((n-1)/n)log2.

The identity sup_{delta>=0}(sqrt(a_n delta)-tau delta)=a_n/(4tau)
then proves the uniform penalty approximation. At a minimizing B_tau,
comparison against a flat optimum gives tau delta<=sqrt(a_n delta),
and hence delta<=a_n/tau^2. These arguments do not assume convergence
or a near-minimizer regularity property.

For the finite-temperature objective, the normalization yields
partial_e^2 f<=beta/n^2, whereas the penalty contributes
-2tau/d=-4tau/[n(n-1)]. Thus its separate concavity threshold is exactly
tau>=beta(n-1)/(4n). Successive coordinate endpoint replacements prove
the claimed equality of the soft and discrete pressure minima. This is
not convexity of the entire objective and licenses no minimax swap.

For the outer Laplace integral, the affine box
(1-epsilon)B_*+epsilon[-1,1]^d has uniform-cube volume epsilon^d.
The strong-variance inequality and Markov give a subset of at least half
that volume on which q(U)<=2sqrt(log2/3). The elementary bound
|B_epsilon^2-B_*^2|<=4epsilon and convexity of Q yield the claimed
objective increment. Integrating only over this subset gives exactly
epsilon(K+4tau)+log(1/epsilon)/alpha+log2/(alpha d). This is uniform
in n>=2 after optimizing epsilon and supports the stated ordered-limit
reduction, but does not itself provide the fixed-parameter limit.

Finally, for C=(A tensor J_k)/sqrt(k), the hollow block magnetization
polynomial reaches its maximum on magnetizations +/-k. Therefore
Q(C)=k^(3/2)Q(A), while its average squared coefficient is
(n-1)/(kn-1). The copying construction consequently retains the cap but
loses the variance occupancy, exactly as claimed.

The sole wording issue reported to the author was to define explicitly
the signed energy q_B(x) used in Section 3, distinct from q_n(B)=Q(B)/n^1.5.
