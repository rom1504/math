# Independent audit of recursive codebook counts

2026-09-07. The mathematical claims in
`flatify_construct_2026_09_07_recursive_codebook_counts.md` pass independent
reconstruction. This audit does not assert a global rank-two cap theorem.

For parent input (xL,xR), the two sparse child words have complementary
supports J,k-J. At fixed J there are binom(k,J)2^k parent inputs, while
the independent signed-support orbit sizes multiply to
binom(k,J)^2 2^k. Their quotient is exactly 1/binom(k,J).
Thus the general histogram overlap formula has the stated normalization.
For fully-flat parent output, both child norms force J=k/2; Walsh
inversion identifies each successful sparse child word with exactly one
1-plateaued Boolean word. This yields the exact squared-count identity.

The cited primary [Potapov Corollary 1](https://arxiv.org/html/2303.16547)
was read directly: its degree bound for 1-plateaued functions is
(d+1)/2. Counting Boolean algebraic normal forms under that degree bound
gives the required exponent k/2+O(k/sqrt(d)). Subtracting log2 binom(k,k/2)
in the overlap formula gives an o(2k) logarithmic upper bound on expected
fully-flat count. Fresh column signs keep every fixed flat tuple's edge
contributions independent symmetric signs of magnitude sqrt(2) in
H_cross/sqrt(N). The stated annealed union bound and rate follow.

For the modular lemma, a sum of an integer word over an r-dimensional
coordinate face is 2^-(d-r) times an integer signed sum of its Walsh
coefficients. If every coefficient is divisible by 2^a and
r>=d-a+1, that face sum is even. In particular the origin-based face
parities giving all Mobius coefficients above degree d-a vanish.
For fixed support, differences of negative-part indicators have Walsh
divisibility 2^(a-1); modulo two they therefore lie in a fixed affine
Reed--Muller space of degree at most d-a+1. Combining its cardinality
bound with the trivial 2^J bound is legitimate, as is intersecting the
support space with weight J via the minimum of its two cardinalities.

For p<=1/2 the resulting per-child-order exponent is
2 min(h2(p),1/2)+1/2+p-h2(p). Below h2(p)=1/2 it is increasing;
above that point it is convex, with the larger endpoint at the transition.
Its maximum is 1+p0, so dividing by parent order 2k gives (1+p0)/2.
All asymptotic Reed--Muller dimensions used here are k/2+o(k), since the
degree cutoff differs from d/2 by only a bounded amount. This is an
upper count for the integral-spectrum sector, not a count for approximate
integrality or all output profiles.
