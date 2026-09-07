# Independent audit: complete fixed-temperature seed-pressure universality

2026-09-07. **PASS.** This audits the complete assembled argument in
`principle_synthesis_2026_09_07_complete_seed_pressure_universality.md`,
including its arbitrary, non-uniformly-integrable row spectra. The descending
matching ingredient was developed jointly; the present audit independently
checks its interaction with the category and conditional-permutation steps.

## Statement and scope

For fixed finite C,K,t, prescribed real row multisets of length m-1 and
individual squared norm at most Cm, and independent uniform row permutations,
write

    Z_S = E exp(-t sum_{i<j}(u_ij-S_ij u_ji)^2).

Uniformly for all such row arrays and Q(S),Q(T) <= K m^(3/2), the proof gives
`|log Z_S-log Z_T|=o(m^2)`. No finite-type, symmetry, minimum atom mass, or
uniform-integrability assumption remains. The per-row energy hypothesis is
essential to the stated matching proof.

## Exact checks

1. The directed-port to edge-difference map has norm at most sqrt(2).
   Coupling the same labeled-copy permutations therefore gives log-partition
   change at most `4t sqrt(C eta) m^2` for squared perturbation at most
   eta m^2. Deleting only selected copies of an equal value still has exactly
   the uniform modified-multiset image law.

2. Descending bins make previously used degree at scale a at most Cm/a^2.
   At a maximal compatible matching, choose one available sign at each
   residual vertex. All unused induced edges are negative after switching.
   The used edges can increase energy by at most twice their number, giving
   `binom(u,2)-u Cm/a^2 <= Q(S)`. Principal cap monotonicity justifies using
   this partial spin. Thus `u <= 2Cm/a^2+1+sqrt(2Q(S))`.
   Summing full row-energy bounds over residual vertices and then the
   geometric bins yields exactly the canonical bound (7). It is loose but
   sufficient and does not assume signs are balanced.

3. Retained matched endpoints belong to the same relative bin and have
   compatible signs. Their squared difference is at most delta^2 times
   half their combined squared magnitude. Hence (9), including its factor
   one-half, holds even when individual values inside a category are later
   freely permuted.

4. The matching prescribes category port sets, not individual values.
   The event probability is exactly the reciprocal category multinomial,
   including one residual category. With reference masses
   `q_b,+ = q_b,- = delta/(8 a_b^2)`, their sum is at most 1/(2V^2).
   Cross-entropy and the row energy bound imply (12); the elementary
   multinomial bound gives (13) with no error depending on the number of
   categories or tiny masses. Arbitrary real values require no extra
   category cost, since within-category ordering remains free.

5. On the category event, replacing retained heavy entries by zero leaves
   exactly the original light multiset together with all the original
   heavy copies turned into zeros. In particular this bounded comparison
   profile is independent of S, despite the preceding unmatched deletion.
   Tagging the required number of zero copies and successively swapping
   each into its prescribed position sends a uniform labeled permutation
   to the uniform conditional labeled permutation. Forgetting labels gives
   precisely the conditional light law. Earlier fixed tags are untouched.
   At most H swaps are used, with H <= Cm^2/V^2; each swap changes at most
   two bounded ports. The canonical 16tL^2 H bound is conservative and
   valid. This is a comparison of normalized conditional expectations,
   so there is no missing probability factor in (15).

6. Each prescribed heavy edge has two heavy endpoints. After collapsing
   them it has kernel one; every other edge is light-light. Thus the
   heavy defect and the conditional bounded partition multiply exactly
   as used in (16). Conversely deleting all heavy entries increases every
   kernel when V>2L, including heavy-light edges. Combining this with the
   deletion perturbation proves both sides of (17).

## Limits and bounded-profile dependency

For a desired accuracy choose delta small, separation ratio R large,
the finite annulus count J large, and then initial radius L0 large.
Only after all these constants are fixed let m tend to infinity. One of
the J disjoint open annuli has energy at most Cm^2/J. Shared endpoints
cause no issue: they are retained and qualify as light or heavy for the
chosen annulus. The chosen index depends only on the row array, not S.
There are finitely many possible bounded radii, so uniformity follows
from the uniform bounded-profile theorem. No diagonal extraction depending
on the seed or hidden growing cutoff is required.

I had previously read the finite-type dependency in full. It permits
row-specific supports in a common finite alphabet. Quantization and
merging rare atoms into a row's most frequent atom retain membership in
that common alphabet. A zero mass is handled by restricting the local
support, while the positive Gaussian kernel dual is smooth on the full
simplex. Thus absent atoms and all-zero rows do not break the dependency.

## Actual-weave application and remaining distinction

Normalized Hadamard spectra have row energy m for every physical spin.
When there is a physical diagonal port, first condition on its assigned
coordinate; the remaining ports still have independent uniform row
permutations and energy at most m. Their profile law is common to the
seeds. Any seed-independent diagonal weight can likewise be incorporated
in this common measure. Uniform multiplicative partition comparison
therefore survives integration over arbitrary common profile-array laws,
including correlated profiles before the final independent permutations.

This removes the non-UI coherent-spectrum escape from the *fixed-temperature
common-law raw-kernel pressure* argument. Exact cycle moments and finite
optimizer-dependent responses remain possible, but cannot by themselves
give a leading m^2 pressure separation here. The theorem does not compare
actual maxima, exclude seed-dependent profile distributions, permit growing
temperature without new estimates, or forbid a genuinely correlated port
placement. It does not establish convergence or the favorable flatification
inequality for selectable optimal children.
