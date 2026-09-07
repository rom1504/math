# Independent audit: stratified marked-selector bulk

2026-09-07. **PASS**, including sharpened clipping and bounded-ratio bands.
**Same-law balanced-face inheritance is separate: L24 remains unproved
by unsigned automorphisms alone; Sylvester L32 avoids that issue.**
Canonical: `principle_invent_2026_09_07_stratified_marked_selector.md`.
This audit does not claim the all-profile cap inequality or convergence.

The exact dephasing identity leaves one marked child and independent ordinary
children: condition on the marked child and absorb its diagonal signs into
each ordinary child's fresh signed input law. The unsigned root permutation
is independent of these displayed children. Its size-L groups therefore
admit independent exactly-k selectors without imposing hidden child
conditioning. Compression of the selector onto group constants is exactly
p times the identity. The relative uniform-selector probability is exactly
binom(L,k)^M/binom(LM,kM), giving the stated entropy fee.

Group words have selector weight binom(L,k)^(-M); type counting gives
H(mu) minus log binom(L,k), not a separate entropy payment in each child.
The child arrays have normalization 1/sqrt(k). Only the marked child has
a raw source, whose mean is zero; the ordinary sources are symmetrized.
The proved marked-spine and ordinary bounds then give Gamma with factor
1/L. All alphabets and depths are fixed before M grows.

The rare-information lemma is independently proved in
`principle_construct_2026_09_07_rare_source_rd_audit.md`. For a packet S,
the nonzero coarse amplitudes are -2 C_j(S)/sqrt(k), with uniformly O(r)
perturbations. Conditional selector entropy is at most log binom(L,k),
so the rare packet entropy is paid once. Its probability lies between r
and kr. No packet-dependent temperature or minimax exchange is used.

The sharper clipping inequality is exact. When the DC term is unclipped,
all terms are unclipped and Parseval gives Ls/(L-1) >= 1+s/(L-1).
Otherwise the DC contributes one. The remaining squared coefficients sum
to s(L-s), and each is at most min(s,L-s)^2. Their clipped sum is at least

    s(L-s)/max(L-1,min(s,L-s)^2) >= s/(L-1),

because both denominator terms are at most (L-1)(L-s). This holds for
every nonempty proper packet, without any support classification.
Consequently tau=k/[4(L-1)] gives the uniform negative coefficient
-k r log(1/r)/[L(L-1)]. Dividing by t times 4r(1-r) gives the cap
coefficient (1-1/L)/(2 sqrt(p)). At k=L-1 this is sqrt(p)/2.

For each fixed K and reference r, all active densities in [r,Kr] obey
the same expansion at the SAME temperature. The uniform O_(L,K)(r)
entropy remainder is o(r_i log(1/r)); constant rows have exactly zero
variance and zero leading exponent. Their polarities cost only exp(O(m)).
Choose child approximation tolerances smaller than the active-row strict
margin and then fix their depths. The macroscopic total-variance union
has only exp(O(m log m)) profile choices. At small total variance the
deterministic operator bound supplies an additive O(delta N^(3/2))
estimate. Sending order first and then delta to zero justifies the
uniform additive-o formulation, not a microscopic relative estimate.

Balanced-column repair is uniform in the spin word. Group independence
gives the needed non-DC excess bound even after stratification. The old
balanced-face certificate is retained only after paying its explicit
selector entropy fee; its separate directed-interval numerical record
is the director's responsibility. Subject to that recorded certificate,
the same ensemble has the two asserted strict-subhalf sectors. Arbitrary
mixtures of balanced and near-constant rows, intermediate means, and
unboundedly many separated minority scales remain unproved.

## Subsequent inheritance caveat

The explicit unsigned Paley12/tensor24 automorphisms were checked exactly,
but they PERMUTE output rows. Averaging a tuple policy under this group
mixes distinct output laws. Invariance of the sum of envelope values
under a permutation before mixing does not prove invariance after mixing.
The envelope is convex, so the needed inequality is not automatic.
Consequently the L24 balanced-face conclusion requires an additional
entropy-compensation proof and is not certified by this audit. For
Sylvester orders, translations act by output SIGNS ONLY and fix the DC
variable exactly. The ordinary child laws are already symmetrized, while
the raw marked DC law remains unchanged. Thus translation averaging
preserves every relevant output law and makes every input marginal the
full balanced source. No preliminary global reversal is needed; such a
reversal could itself mix the raw marked law and must not be assumed
harmless. This proves the safe Sylvester inheritance step.
