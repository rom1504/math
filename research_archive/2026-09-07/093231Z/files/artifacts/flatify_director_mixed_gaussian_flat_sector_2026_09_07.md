# A simultaneous Gaussian/flat profile theorem for actual rank-two signs

2026-09-07. Proved by the director and independently reconstructed by the
adversarial researcher. This is a quantitative sector theorem, not an
all-spin cap improvement or a selected-child recovery theorem.

Use the actual rank-two cross construction with k=2m, N=mk and a fixed
order-k Hadamard frame, independently signed/permuted in each fibre.
For a physical spin x let rho_i be the empirical magnitudes of its
normalized transformed row. Before the random permutations, assign each
row either label G (Gaussian) or label B (fully flat). Write b for the
fraction of B rows and g=1-b, and set

    delta(x)^2 = (1/m) sum_i W2(rho_i,nu_i)^2,
    nu_i = law(|standard Gaussian|) for G, and delta_1 for B.

The assignment may depend on the whole physical spin x; it must not depend
on the random frame realization. No claim that B rows are exactly bent is
needed. Exact flat spectra are simply the zero-distance case.

Put cG=sqrt(15)/8 and kappa=5/(16sqrt(15)). Choose epsilon such that the
empirical law of k independent half-normal samples is within W2 distance
epsilon of the half-normal law with probability at least1/2. There exists
ONE realization of the cross signing such that, for EVERY physical spin,

    |H_cross(x)| / N^(3/2)
      <= cG-kappa*b(x)^2
         +cG*g(x)*(2epsilon+epsilon^2)
         +(delta(x)+epsilon)*(2+epsilon)/2.               (1)

One may use epsilon=O(k^(-1/4)sqrt(log k)). Thus any sector satisfying
delta(x)=o(1), with arbitrary mixtures of the two row classes, has uniform
cap at most cG+o(1). Completing the fibres with actual low-cap signs adds
O(N^(5/4)). No Boolean row-count bound, independent Gaussian replacement
of the OUTPUT signing, or separately paid scalar channel is used.

## Proof

The cross micro-operator K has norm1. Couple each G row to independent
Gaussians conditioned on its empirical-magnitude W2 event, exactly as in
`flatify_adversary_2026_09_07_gaussian_profile_sector.md`. Couple each B row
to independent fair signs of magnitude1, using the same random signs and
ordering as its physical word. With u the actual feature vector and v the
comparison vector, sorted matching gives

    ||u-v|| <= (delta+epsilon)sqrt(N),
    ||u||=sqrt(N), ||v||<=(1+epsilon)sqrt(N).

The quadratic error is at most eN, e=(delta+epsilon)(2+epsilon). The
Gaussian coordinates alone have squared norm at most g(1+epsilon)^2N.
If there are s=gm Gaussian rows, their joint conditioning costs at most2^s.

Use t=2sqrt(15), a=16, lambda=(a-1)/2=15/2. For the unconditioned
comparison variables insert exp(-lambda times Gaussian squared norm).
Independent edge-coordinate groups give these EXACT factors:

    GG edge: (a^2-4t^2)^(-1) = 16^(-1);
    GB edge: a^(-1) exp(4t^2/a);
    BB edge: cosh(2sqrt(2)t) = cosh(4sqrt(30));
    two unused coordinates of one G row: a^(-1).

For GB, the rotated flat pair has norm sqrt(2), so the ordinary Gaussian
linear moment gives the displayed exponent. For BB the Boolean H2
identity gives precisely two energy values, not independently bounded
channels. There are s(s-1)/2 GG edges, s(m-s) GB edges and
(m-s)(m-s-1)/2 BB edges. Let A=log16, L=logcosh(4sqrt(30)). The leading
logarithmic moment INCLUDING the norm compensation, divided by N, is

    F(b)=lambda*g-g^2*A/4+g*b*(-A+15)/2+b^2*L/4
        =F(0)+(b^2/4)*(A-30+L),
    F(0)=15/2-log2.

The linear term cancels exactly. Also

    A-30+L < 3-30+22 = -5,

using log16<3, logcosh(z)<z and sqrt(30)<11/2. Since
log2+F(0)=15/2=2t*cG, the leading exponent in the two-sided full-cube union
bound is nonpositive at (1). Gaussian norm enlargement contributes
lambda*g*(2epsilon+epsilon^2)N; the coupling contributes t*eN. Both are
explicitly paid in (1).

The exact finite correction, including conditioning and both polarities,
is

    log2+s log2-s*A/2-(m-s)*L/2 <= (1-m)log2.

Here L>2log2. Union over at most2^N physical spins therefore leaves failure
probability at most2^(1-m)<1 for m>=2. Each spin can have its own labels
and coupling: the required probability estimate is fixed-spin before the
union bound. This proves simultaneous existence. The elementary empirical
transport rate and actual within-fibre completion are detailed in the
Gaussian-sector artifact.

## What was learned, and what is still missing

Joint heterogeneous interactions do not spoil the Gaussian-sector bound
when the other rows are fully flat or close to it. In fact they produce a
strict quadratic improvement in the flat-row fraction. This is stronger
than separately proving two homogeneous sectors or paying their cross
interaction by a norm triangle bound.

The theorem does not cover every profile: exact sparse spectra, intermediate
non-Gaussian spectra and certain mixtures still require other arguments.
The exact-flat-spectrum theorem controls the homogeneous sparse class by
different arithmetic. Combining those results does not automatically control
their mixed interactions with all remaining profiles. Nor does this random
construction preserve its outer seed. Convergence remains open.

The exploratory optimizer in
`computations/flatify_director_gaussian_bent_mixture_2026_09_07.py` suggested
the cancellation, but all constants in (1) follow from the analytic proof
above. Its floating-point outputs are not used as certificates.
