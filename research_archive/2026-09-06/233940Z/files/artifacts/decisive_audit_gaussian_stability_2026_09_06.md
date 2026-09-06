# Independent audit: Gaussian deformation and deficit confinement

Date: 2026-09-06. Status: PASS, at the exact scope stated in
`decisive_independent_gaussian_stability_2026_09_06.md`.

I independently reconstructed both proofs rather than relying on their
reported verdict. Neither establishes convergence, and neither asserts
exponentially many exact ground states.

## 1. Deformation normalization and error

Let `a=(1+t²)^(-1/2)` and `sigma²=1-a²`. Independent signs with means
`a A_e` have centered variance `sigma²`; their centered absolute third
moments are at most 8. The Gaussian replacement third moment is less than 2.
For the displayed soft maximum, an edge third derivative is a third
centered moment of a variable in `{−1,1}`, multiplied by
`beta²/n^(5/2)`, and hence has absolute value at most
`8 beta²/n^(5/2)`. Taylor replacement therefore costs at most
`(20/3) beta²/sqrt(n)` after summing fewer than `n²/2` edges.

Only ONE soft-maximum error is needed. The sign-side soft maximum is
already at least the sign-side maximum, hence at least `m_n`; the Gaussian
side maximum is at least its soft maximum minus
`(n+1) log(2)/(beta n)`. With `beta=n^(1/6)` the sum of constants is less
than 10. Finally `a A+sigma G=a(A+tG)`, which gives precisely

`F_A(t) >= sqrt(1+t²) (m_n-10 n^(-1/6))`.

All comparison estimates are uniform in the deterministic signing. This
uses actual global minimality through the randomized *sign* matrix, not
an unproved minimization/expectation exchange.

## 2. Deficit-shell upper bound

For projective Hamming distance `r`, the Gaussian quadratic increment
has variance exactly `4r(n-r)`. The extra factor 2 in the shell counting
is harmless and covers both representatives of a projective spin.
The maximum over a fixed shell is Gaussian Lipschitz with constant at
most `n`, since every coefficient-difference vector has Euclidean norm
at most `n`. Gaussian concentration and a union bound over at most `n+1`
shells give the stated `n sqrt(2 log(n+1))` error before multiplication by
the fixed parameter `t`; its normalized contribution vanishes.

The center maximum contributes asymptotically at most `t sqrt(eta)`.
The shell term is bounded by

`2t sqrt(eta) + sup_(0<=rho<=1/2)
 rho[-kappa+sqrt(8)t sqrt(log(e/rho))]`.

For a positive value, `L=log(e/rho)>L0=kappa²/(8t²)>=2`.
Dropping the negative term bounds it by
`sqrt(8)t exp(1-L) sqrt(L)`, decreasing on this interval; its value at
`L0` is exactly `e kappa exp(-kappa²/(8t²))`. Combining with Section 1
proves the displayed entropy-rate inequality with the claimed constant 3.

## 3. Quantifiers and limitations

The result applies to any asymptotically minimizing sequence and any
subsequence on which its normalized optimum and center entropy rate
converge. A positive asymptotic lower bound for `m_n` is enough; the
elementary bipartition/Khintchine argument gives `1/4+o(1)` with the
quadratic normalization used here.

The confinement assumption concerns EVERY signed configuration and
its deficit from the global absolute cap. It is much stronger than saying
that exact grounds have a small covering number. No conclusion about
ground-state cardinality, exact minimizer classification, or cross-order
realizability follows without an additional argument. This distinction is
correctly retained in the source artifact.
