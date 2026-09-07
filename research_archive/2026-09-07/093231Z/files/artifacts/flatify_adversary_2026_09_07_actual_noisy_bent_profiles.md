# Actual growing-order noisy-bent row profiles

2026-09-07. A proved realizable family and an entropy lower bound. It is
not a classification of high-entropy rows or an upper bound on all type
counts. The numerical transport experiment below is not a certificate.

Let k=4^r and F be the Walsh Hadamard matrix of order k. Fix a bent sign
row b, for example the sign of the bilinear form on two r-bit halves.
Then c=F^T b/sqrt(k) is another sign row. For fixed epsilon in (0,1/2],
put alpha=1-2epsilon and take

    x_i=b_i eta_i,
    P(eta_i=-1)=epsilon,

independently in i. Write z=F^T x/sqrt(k). After multiplying coordinate j
by c_j, its mean is alpha and its centered part is a sum of k independent
bounded variables with coefficients ±1/sqrt(k). Its variance is
1-alpha^2. Distinct centered coordinates have covariance zero, by
Hadamard orthogonality.

The bounded-coefficient one- and two-dimensional Lindeberg central limit
argument is uniform over coordinates and distinct coordinate pairs. Thus,
for each bounded Lipschitz test function f, the mean of
k^(-1) sum_j f(|z_j|) tends to the corresponding folded-normal expectation
and its variance tends to zero. It follows that the empirical magnitude
law converges in probability to

    rho_epsilon = law(|alpha+sqrt(1-alpha^2)G|).

Its second moment is one, as is the empirical second moment for every
physical x. Uniform bounded fourth moments (expand the independent centered
sum) control tails, so convergence also holds in W_2 in probability.

## Physical entropy and actual types

The product-noise typical set has cardinality

    exp[k h(epsilon)+o(k)],
    h(epsilon)=-epsilon log epsilon-(1-epsilon)log(1-epsilon).

Intersecting with a shrinking W_2 neighborhood of rho_epsilon retains
probability tending to one. Restrict at the same time to a shrinking
Hamming-weight window around epsilon k. Each remaining word then has
probability exp[-k h(epsilon)+o(k)], uniformly, so the intersection still
contains exp[k h(epsilon)-o(k)] distinct PHYSICAL Boolean rows.

There are only exp(O(k^(2/3))) possible Hadamard magnitude types, by the
square-partition count in the audited typed certificate. Pigeonholing gives
at least one realizable type rho_k with

    W_2(rho_k,rho_epsilon)->0,
    liminf k^(-1) log c_k(rho_k) >= h(epsilon).

This lower entropy bound is genuine and can be used to challenge a proposed
uniform row-profile theorem. For epsilon=1/2 it attains log2, the universal
upper bound on physical entropy. For other epsilon, c_k(rho_k) may be
larger because rows outside this one noisy-bent family can have the same
type. No upper bound h(epsilon) on the full c_k is inferred.

Different epsilon values can coexist at the same k in different fibres;
thus heterogeneous tuples of these realizable families are available.
No transport-continuity theorem at unbounded exponential costs is being
inferred merely from the W_2 convergence above.

## Numerical diagnostic, with its precise limitation

`computations/flatify_adversary_2026_09_07_noisy_bent_transport.py` evaluates
the homogeneous four-coordinate transport dual for 12- and 20-point
Gauss-Hermite discretizations of the folded normal laws, at t=1,2,4,8,16.
The rate inserted into the displayed expression is the proved LOWER bound
h(epsilon), not an assumed value of the full type entropy.

At 20 points the minimum sampled expressions for alpha=0,.2,.4,.6,.8,.95
were approximately .484123,.482874,.478991,.465405,.429306,.333582.
No counterexample was found by this diagnostic. In particular this does
NOT prove that every type converging to a noncentral folded Gaussian is
controlled; its full physical entropy and a rigorous transport passage
would still be needed. All supports, weights, potentials, optimizer
residuals, and outputs are saved in the same-stem result JSON.

For the central Gaussian law the transport is analytically solvable, and
the separate actual Gaussian/W_2 sector theorem is proved in
`flatify_adversary_2026_09_07_gaussian_profile_sector.md`. That theorem uses
the full 2^N spin count and does not rely on this numerical experiment.
