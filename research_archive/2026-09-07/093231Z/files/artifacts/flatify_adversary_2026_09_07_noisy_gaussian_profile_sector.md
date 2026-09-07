# Uniform heterogeneous noisy-Gaussian profile sector

2026-09-07. A proved extension of the Gaussian/flat sector. It does not
claim to cover all physical Hadamard row types.

Let R be a fair sign and G standard Gaussian, independently, and define

    X_alpha=alpha R+sqrt(1-alpha^2)G,  0<=alpha<=1,
    rho_alpha=law(|X_alpha|).

For the actual rank-two cross-sign ensemble, assign to every physical
spin row any alpha_i before the random signed column permutations. Let
delta^2 be the average squared W_2 distance of its row magnitude law to
rho_(alpha_i). Then there are actual outputs controlling ALL such spins
simultaneously by

    Q_sector/N^(3/2) <= c_NG+delta+o(1),
    c_NG=sqrt(15)/8+log(16/15)/(16sqrt(15))
        =.4851644041824049... .

The o(1) is uniform in every choice of alpha_i and bounded delta. No
physical-entropy assignment h(epsilon) is used: the union counts all 2^N
spins. Low-cap within-fibre completions add O(N^(5/4)).

## Exact tilted joint kernel

Use t=2sqrt(15), lambda=15/2. Take independent coordinate sources
X_alpha,X_alpha at one endpoint and X_beta,X_beta at the other. Include
the tilt exp(-lambda sum of their four squares). Conditional on their
Rademacher means, integrate the Gaussian parts. If p=alpha beta and
D=16-15p^2, the exact resulting edge expectation is

    M(p)=D^(-1) exp(-30p^2/D) cosh(4sqrt(30)p/D).           (1)

Here is a direct normalization check. Write O=H_2/sqrt(2), so the edge
quadratic is 2X^T OY. The two Gaussian variances are sigma^2=1-alpha^2 and
tau^2=1-beta^2. The precision determinant per coordinate pair is

    (1+2lambda sigma^2)(1+2lambda tau^2)-4t^2sigma^2tau^2
       =16-15alpha^2beta^2=D.

After integration the effective diagonal coefficients on the two mean
vectors are -lambda beta^2/D and -lambda alpha^2/D, and the cross
coefficient is t O/D. Their mean vectors have squared norms 2alpha^2 and
2beta^2. Finally a two-bit sign pair satisfies R^T O R'=±sqrt(2), giving
exactly (1), including its factor four in the cosh argument.

The unused single-coordinate factor is

    J(alpha)=(16-15alpha^2)^(-1/2)
             exp[-(15/2)alpha^2/(16-15alpha^2)] <=1/4.     (2)

Indeed, the derivative of log J(alpha)^2 with respect to r=alpha^2 is
-225r/(16-15r)^2. Thus it decreases from J(0)^2=1/16.

## Rational uniform certificate for M

The executable certificate proves M(p)<1/15 throughout [0,1]. Divide this
interval into 500 rational cells [l,u], set Dmin=16-15u^2 and
Dmax=16-15l^2, and use the verified rational upper square root
R=2738612787525831/500000000000000>sqrt(30). Monotonicity gives

    M(p) <= [exp(S-A)+exp(-S-A)]/(2Dmin),
    A=30l^2/Dmax,  S=4Ru/Dmin.

Every one of these 500 upper bounds is strictly less than 1/15 by exact
rational arithmetic. The maximum leaf upper is

    11163527888789114680833563530151267429109
    /171315200000000000000000000000000000000000.

The exponential routine range-reduces by 128, uses the degree-18 Taylor
polynomial and the absolute geometric remainder
|z|^19/(19! (1-|z|/20)), rounds upward to denominator 10^40, then squares
upward seven times. All reduced |z|<1. Acceptance uses no floating
arithmetic. All rational cell bounds are saved.

- `computations/flatify_adversary_2026_09_07_noisy_gaussian_kernel_certificate.py`
- `computations/results/flatify_adversary_2026_09_07_noisy_gaussian_kernel_certificate.json`

## Coupling and all-spin estimate

Use the sorted-coordinate coupling from
`flatify_adversary_2026_09_07_gaussian_profile_sector.md`, now sampling each
row from its own X_(alpha_i) law conditioned on empirical magnitude W_2
distance at most epsilon. These laws are uniformly subgaussian with
variance one, so the same truncation/CDF proof gives probability at least
1/2 at epsilon_k=O(k^(-1/4)sqrt(log k)), uniformly in alpha.

The actual signed-permutation features u and comparison features v satisfy

    |u^T K u-v^T K v| <= eN,
    e=(delta+epsilon)(2+epsilon),
    ||v||^2 <= (1+epsilon)^2 N.

Unconditioned comparison coordinates are independent, so (1)--(2) give
moment at most 15^[-m(m-1)/2]16^(-m). Conditioning costs at most 2^m.
Markov and the 2^N spins plus two polarities imply a bad-event bound whose
leading exponent per N is

    log2+lambda-(1/4)log15-2tq
      +t e+lambda(2epsilon+epsilon^2).

The remaining finite logarithmic correction is

    m[(1/2)log15-log16+log2]+log2 <=(1-m)log2.

Therefore the finite choice

    q=c_NG+e/2+[lambda/(2t)](2epsilon+epsilon^2)

has bad-sector probability at most 2^(1-m)<1. This proves the asserted
uniform sector theorem and retains the exact loop factors.

## Explicit limitation of the reference family

This kernel bound does NOT extend unchanged to all normalized Rademacher
sums. Already X=(R_1+R_2)/sqrt(2), with magnitude law supported on 0,sqrt(2),
has tilted self-edge expectation approximately .0817011660>1/15. That
source is handled by the separate plateaued-sector bound, but it refutes
blindly applying (1) or its uniform constant to multiple Rademacher means.
