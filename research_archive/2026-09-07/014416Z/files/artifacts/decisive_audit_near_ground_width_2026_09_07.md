# Independent audit: deficit-conditioned Gaussian width

Date: 2026-09-07. Verdict: PASS for the exact statement in
`decisive_director_near_ground_width_susceptibility_2026_09_07.md`.
This is not a convergence theorem.

I independently checked the deformation input in
`decisive_independent_gaussian_stability_2026_09_06.md`, rather than
assuming its earlier audit. With the normalized soft maximum
`(beta n)^(-1) log sum exp(beta H/sqrt(n))`, an edge third derivative
is bounded by `8 beta^2/n^(5/2)`. The centered biased sign and matched
Gaussian have total absolute third moment below 10. Fewer than `n^2/2`
Lindeberg steps and the Taylor factor `1/6` give at most
`(20/3) beta^2/sqrt(n)`. The soft maximum error is at most
`2 log(2)/beta`. At `beta=n^(1/6)`, the stated `10 n^(-1/6)` is safe.
Dividing by the bias `a=(1+t^2)^(-1/2)` gives the claimed factor
`sqrt(1+t^2)` including its error. This input genuinely uses minimum
over actual signings, not convex optimality of a weighted matrix.

For fixed bin width d, there are finitely many deficit bins, uniformly
in n along the selected order subsequence. Every bin Gaussian maximum
is `sqrt(binom(n,2))`-Lipschitz. The expectation/max interchange costs
at most `t sqrt(n(n-1) log J)=o(n^(3/2))` at fixed J. Thus finite maxima
commute with the necessary limsup in the valid direction, giving

`c(sqrt(1+t^2)-1) <= max_j[-j d+t W((j+1)d)]`.

Suppose the desired limsup were strictly smaller than `2c`. Choose
`K^2<2c` and e0 with `W(e)<=K sqrt(e)` for ALL `0<e<=e0`.
Choose t small enough that bins beginning beyond e0/2 have negative
value using the global width bound, then take `d<e0/2`. The remaining
bin maximum is at most `K^2 t^2/4+d`. Let d tend to zero AFTER the order
limit and then t tend to zero. This forces `c/2<=K^2/4`, a contradiction.
No uniformity in a moving finite-n deficit window is used or concluded.

The result is not a cardinality tautology. Fix one Boolean center z and
consider its projective Hamming ball of radius rho n, with one energy
orientation. At distance r the Gaussian increment has variance
`4r(n-r)`. A Gaussian maximum bound over the ball gives normalized width
`O(rho sqrt(log(e/rho)))`, after n tends to infinity. This is much
smaller than `sqrt(rho)`, despite exponentially many configurations at
each fixed rho>0. Thus a landscape whose e-near-ground states all lie
in finitely many `O(e)n`-radius balls would not automatically satisfy
the theorem. This illustration concerns geometry; it is not asserted
to construct such a landscape for actual minimizing signings.

The theorem permits `W(0)>0`, when its conclusion is immediate. It
does not imply exponentially many exact ground states, an overlap law,
transversality to deterministic completion matrices, or an all-order
comparison.
