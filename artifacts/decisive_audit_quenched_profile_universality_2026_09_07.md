# Independent audit: quenched finite-profile universality

Date: 2026-09-07. Full proof audit PASS for
`decisive_director_quenched_finite_profile_universality_2026_09_07.md`.

The hypothesis `Q(S_n)=o(n^2)` supplies `beta(S_n)=o(n^2)` by the
elementary polarization inequality. Indicators of arbitrary vertex
classes are allowed bounded vectors, so the class-pair sign frequencies
are asymptotically balanced uniformly over class assignments, including
assignments depending on `S_n`. A within-class ordered sum counts edges
twice and loses only the hollow diagonal; the stated coefficient `1/4`
correctly includes the unordered-edge factor and the fair-sign factor.

The lower bound does NOT need every vertex to have balanced neighborhoods
inside every class. For the independently drawn edge colors, the sum over
all vertices of each class-neighborhood discrepancy is bounded by
`||S_n 1_(V_b)||_1<=beta(S_n)`. Summing finitely many classes and colors
gives total expected row-count error `o(n^2)`. Proportion/profile errors
and independent-edge fluctuation errors have the same negligible order.
This is exactly the quantity required by endpoint recoloring: a few bad
rows may cost `O(n)` each without invalidating the total repair budget.

Both fixed positive kernels have bounded logarithms. Therefore endpoint
repair changes log weight by `o(n^2)` and has `exp(o(n^2))` preimage
multiplicity. The independent-edge typical-set information density uses
only the fixed nonzero entries of the chosen couplings; zero entries are
never drawn. This completes the lower bound at the stated quantifiers.

For the centrally symmetric homogeneous case, represent a feasible pair
of couplings as a law of `(e,X,Y)` with fair `e`. The two transformations
`(e,X,Y)->(-e,bar X,Y)` and `(-e,X,bar Y)` generate a four-element group.
The kernel identities preserve cost, and the fixed unconditioned marginal
`nu=bar nu` preserves feasibility. Averaging the WHOLE group leaves the
conditional marginals of both coordinates equal to
`(nu+bar nu)/2=nu` in each sign sector. It also preserves endpoint-swap
symmetry, because swapping endpoints conjugates the two generators.

Reflecting the negative sector's second coordinate unfolds its cost to
the positive kernel. Its two marginals remain `nu`; if the resulting
coupling is not symmetric, endpoint symmetrization increases entropy and
preserves cost. This justifies the common self-transport upper bound.
For the reverse inequality, simultaneously-reflection-invariant and
swap-invariant positive optimal coupling exists by averaging. Its
one-coordinate reflection is symmetric as required for the negative
sector and attains the same value. Fixed colors of the involution do not
affect any step. The final coefficient is `1/2`, as stated.

This is an exact fixed-alphabet pressure limit for each qualifying outer
sequence. It gives no uniform two-sided approximation for unbounded
spectra: a vanishing number of incidences can carry a nonvanishing share
of the quadratic defect. Consequently it cannot identify the cap of an
actual optimizing sequence or erase seed dependence residing in energetic
spikes. The source explicitly preserves that distinction.
