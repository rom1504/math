# Independent audit and simplification of entrywise-small surgery

Date: 2026-09-07. Verdict: the constructive entrywise surgery and its
minimax covariance consequence are correct. Target contraction removes
the original bias term and the need for an input cap bound.

Let A be a hollow full signing, |G_ij|<=1, and 0<t<=1. Set

    B=(A-tG)/(1+t).

This is feasible as a sign mean. Scalar contraction gives
Q(B)<=Q(A-tG), exactly. Independent sign rounding has edge flip
probability t(1+A_ijG_ij)/[2(1+t)] and centered edge variance at most

    1-[(1-t)/(1+t)]^2=4t/(1+t)^2<=4t.

With a=(n+2)log2, the full-cube Bernstein argument consequently gives
an actual A' with

    Q(A')<=Q(A-tG)+sqrt(4t n(n-1)a)+4a/3.

For t=theta/sqrt(n), the normalized error is

    e_n=2sqrt(theta(1-1/n)a/n)n^(-1/4)+4a/[3n^(3/2)].

The earlier extra theta C/sqrt(n) is unnecessary. The edit-count Markov
argument survives with an improved expectation. No rank, PSD, coherence,
or cap-size assumption on A is needed for this actual sign construction.

At a normalized delta-near-minimizer, uniformity over the entire symmetric
entrywise cube permits minimax. Its support function is the sum of the
absolute off-diagonal coefficients, giving one signed-spin law with

    E[q-sigma h]+(theta/n^2)sum_(i<j)|E[sigma x_i x_j]|
      <=delta+e_n.

The normalization theta/n^2 is correct: the perturbation coefficient is
t=theta/sqrt(n) and energies are divided by n^(3/2). There is no factor
two because energies sum over unordered edges. Diagonal entries are
excluded, so the bound does not constrain E sigma.

If the negative endpoint is gamma below q, opposite-polarity mass is at
most (delta+e_n)/gamma. The stated positive-polarity off-diagonal bound
then follows from triangle inequality and at most n(n-1)/2 entries.
This is a spread constraint on near-ground laws, not an added-row bound
or a proof of balanced cap endpoints.
