# Independent audit: correlated entropy extraction

Date: 2026-09-07. PASS for the full finite statement in
`decisive_independent_correlated_entropy_extraction_2026_09_06.md`.

The gauge and orientation placing an absolute extremizer at the all-one
spin preserve the actual cap. The simultaneous two-sided diagonal
majorant gives a core of size at least `n/2` with normalized operator
bound `K=8K_G C`; no claim is made that the whole parent has this bound.
Choosing `kK<=1/2` makes the core Gaussian covariance positive definite.
The outside coordinates remain present, independently biased, so the
entropy baseline is `n h(delta)`, not just the core entropy.

The Gaussian relative entropy is exactly `-(1/2)logdet Sigma`, since
the covariance diagonal is one. For `|u|<=kK`, integration of
`u/(1+u)` gives `-log(1+u)+u<=u^2/[2(1-kK)]`. The hollow trace cancels
the linear term, and `tr B_I^2=m(m-1)`. Thus the entropy cost is exactly
bounded by `k^2m(m-1)/[4n(1-kK)]`, with the stated factor.

For the biased threshold, orthonormal Hermite covariance is
`r^2+sum_j c_j^2 u^j`, with `c_1^2=4phi(z)^2` and total remaining mass
`1-r^2`. Multiplication by a flat edge sign makes every odd term
positive, independently of that edge sign. Every even term has the SAME
coefficient on every core edge and hence sums against `H_(B_I)(1)`.
Principal monotonicity bounds its absolute value by `Q(B)`, producing
only `beta C(1-r^2)k^2` in log-partition error. This exact cancellation
would be lost under an absolute edgewise remainder estimate.

The resulting extensive coefficient is
`2 beta a k-k^2/[4(1-kK)] >= (3/2) beta a k`, because
`k<=beta a` and `kK<=1/2`. Finally
`m(m-1)/n>=n/4-1/2`. These give the displayed finite inequality, not
merely its asymptotic form. The Gibbs variational inequality uses the
chosen orientation's partition function, which is bounded by the sum
of both orientations, so no factor two in the energy is missing.

The strictly positive pressure correction is uniform over actual
cap-bounded parents for each FIXED positive beta and interior delta.
The proof does not prove a uniform gap as delta approaches its boundary,
or compactness of optimized corrected parameters. Its stated limitation
on the infimum over all product-corrected choices is therefore necessary.
