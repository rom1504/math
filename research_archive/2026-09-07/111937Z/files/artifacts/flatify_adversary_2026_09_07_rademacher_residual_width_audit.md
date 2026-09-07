# Independent audit: Rademacher residual width

2026-09-07. PASS, with the original certificate-only scope retained.

Audited `flatify_independent_2026_09_07_rademacher_residual_width.md` from its definitions, rather than importing the earlier Gaussian-width verdict.

For the deterministic split |P|~2n/3, |Q|~n/3, the greedy difference is exactly 2Z, with Z=sum_i |(B eta)_i|. Therefore E Z/n^(3/2) tends to (2/3)sqrt(2/(3pi)); the normalization is half-width, not absolute cap at one designated polarity. Z<=Q(A) follows because the two energies differ by 2Z and each has absolute value at most Q(A).

The cap-to-fourth-trace estimate is valid: use the i-th row of A as a cube vector to get sum_j |(A²)_ij|<=beta(A), and |(A²)_ij|<=n. Summing gives tr A^4<=n² beta(A). Restriction B=P_P A P_Q contracts the Schatten-four norm, so ||BB^T||F²<=tr A^4.

The asserted two-row correlation estimate requires no multivariate CLT. Gauge the first row to ones. If the agreement and disagreement counts are a,b, the two fields are U+V,U-V, with independent sums U,V. Away from field ties their sign product is sign(|U|-|V|). Scalar Kolmogorov error controls the absolute-value CDF with twice the error; integrating against an independent second marginal changes the comparison probability by at most the sum of these errors (plus atoms at the boundary). When |rho|<=1/2, a,b>=q/4, and field ties have probability O(q^-1/2), including both parity possibilities. Gaussian comparison gives (2/pi)arcsin rho. Thus the uniform bound C(|rho|+q^-1/2) is justified. The sole external input is the standard scalar bounded-summand Berry--Esseen theorem stated explicitly in the source. Its linked publisher page could not be fetched during this audit; this is not a claim to have independently read that paper.

Squaring and summing the correlation estimate yields O(n^1.5) for the squared Frobenius norm of the greedy-sign covariance. The mixed block identity E eta_j s_i=B_ij E|S_q|/q is exact even for even q because fair ties have zero conditional expectation. Its squared Frobenius norm is at most |P|. Consequently ||K||F=O(n^.75).

For an arbitrary rank-r projection, tr(P_X K)<=sqrt(r)||K||F. Hence r=o(sqrt(n)) suffices to make the projection energies o(n) in expectation. Removing a vanishing-probability exceptional event costs at most Q(A) times its probability from E Z, since Z is nonnegative and bounded by Q(A). This proves the claimed uniform near-subspace half-width constant.

Finally, the two child energy differences may be combined without choosing their individual signs: the difference between the plus/plus and minus/minus internal sums is twice the sum of the half-widths. One of the two absolute sums is at least that sum. Both residual norms are asymptotically sqrt(n). Thus the separately paid scalar envelope is at least [1+(4/3)sqrt(2/(3pi))-o(1)] n^(3/2).

No lower bound on the actual bridge energy at those chosen child witnesses follows. In particular, this result does not refute the new localized full-sign bridge operation or unrestricted global flatification.
