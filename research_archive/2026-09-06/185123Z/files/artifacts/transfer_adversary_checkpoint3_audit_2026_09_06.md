# Independent audit handoff, checkpoint 3

Date: 2026-09-06, approximately 18:49 UTC. This records completed
adversarial audits and remaining scope limits. It is not a proof of
convergence or nonconvergence of the original normalized minima.

## Completed new proofs owned by this track

- `transfer_adversary_fixed_retention_random_loss_2026_09_06.md`:
  for each bounded parent-cap coefficient and target `c<2/pi`, a
  sufficiently small positive retention interval has uniform low-cap
  selector probability at most `K/n`. The covariance proof treats actual
  and ambient coloring overlaps separately; either nonempty overlap
  merges two connected Eulerian components and gains exactly `1/n`.
  The checker passed 200 exact ambient identities, 800 subset covariance
  cases, and 4384 connected-overlap exponent checks.
- `transfer_adversary_deterministic_planting_escape_2026_09_06.md`:
  any prescribed child `B_n` can be placed exactly in a parent differing
  only on that principal block, at finite cost
  `Q(A') <= [1+2 K_G p/(1-p)]Q(A)+Q(B_n)`, where `p=n/D`.
  In particular every prescribed bounded-cap sublinear child occurs
  inside an asymptotically minimizing parent. This does not assert that
  an exact minimizing parent already contains it.

The planting operation uses a spectral core solely to locate the old
block. It does not remove the remaining parent vertices. Thus fixed
deletion `epsilon=1/2` gives a loss of order `p`, not `sqrt(p)`.

## Independent audits of other tracks: PASS

### Exponential selector cost

`transfer_director_exponential_selector_cost_2026_09_06.md` was checked
including its primary convex-distance mapping, finite constants,
core reduction, sampling estimate, and relative-entropy consequence.

The active positive-semidefinite quadratic has subgradient norm at
most `2 L sqrt(n)` on the exact binary slice. Proposition 1.7 of
Sambale--Sinulis, arXiv:2010.16289v1, therefore gives the denominator
`144*4=576`. The primary proposition imposes no permutation-invariance
condition on the event; the different sampling-space theorem elsewhere
in that paper is not being substituted for it.

For the ordered sample without replacement, the conditional range
of the Doob martingale increment is exactly
`1-(n-k)/(D-k)=(D-n)/(D-k)<=1`. The bounded-interval exponential
estimate gives `E exp(t(m-E m))<=exp(nt^2/8)`, hence the core-size
tail `exp(-2 epsilon^2 n)`. Conditional on its size, the intersection
is uniform in the core. The cap-gap exponent and final factor three
are correct. Binary relative entropy discards only the nonnegative
term `-(1-gamma)log(1-U(G))`.

The finite checker passed all 139264 tangent-pair cases over the 64
order-four hollow signings. Minor clarity suggestions were sent to the
owner: explicitly call the bounded-operator input a hollow signing;
call its local error estimates mean-square estimates; and, if desired,
choose the small retention constant at most one half to exclude the
full slice automatically. These are not substantive proof gaps.

### Finite exceptional-selector defects

`transfer_seed_exceptional_selector_finite_defects_2026_09_06.md` was
checked throughout. Covariance clipping creates an actual correlation
matrix even without a positive diagonal floor. The two discarded row
sets cost at most `2 delta sqrt(s)/b`, with the stated hollow-energy
normalization. The two-sided cubic calculation cancels odd global
moments, and optimization in `b` gives exactly

`Delta >= [Gamma_+^2/(8 pi c sqrt(S))-sqrt(S)/(2 sqrt(n))]_+`.

The rational limiting defect `841/2213750>1/3000` and the finite
screen `0.5023561171593656...>1/2` check. The replay passed 1120
clipping and 280 two-sided cubic cases, as well as its exact constants.

### High-retention / sparse-retention transition

`transfer_seed_random_retention_transition_2026_09_06.md` was checked
against the actual standalone strict-upper proof, equations (28) and
(34)--(36), not merely its new summary. Input signed-permutation
invariance gives the same exponential failure estimate for every fixed
balanced selector in one ensemble. Averaging that fraction genuinely
selects one deterministic full Hadamard parent.

Conditional on any fibre-count vector, the repair operation has the
product-uniform balanced law, independent of those counts. The exact
hypergeometric variance and cap continuity bound give its stated
finite probability estimate. The same weave has actual Boolean
eigenvectors for both eigenvalues `+m` and `-m`, so its normalized
full hollow cap tends to one half. Consequently the same deterministic
parents have typically strict-subhalf children at retention `31/32`
and typically at least `2/pi-o(1)` children at every vanishing retention.

The replay passed 1994 exact coupling cases and 100 cap-continuity
checks. A harmless `u` versus `nu` symbol typo was sent to the owner.

## Scope and genuinely open work

Two different scope limits are now proved, not conjectured. Small
retention is essential to the universal random-selector loss, because
high retention can be typically good. Random-selector rarity also
does not exclude deterministic transfer: specially planted sparse
selectors can be good even inside asymptotically minimizing parents.

Neither construction transports an arbitrary small minimizing seed
upward at its original coefficient. The original convergence question,
and the independently studied stabilized Hadamard catalyst equality
`R=T`, remain open in this work. The order-31 arithmetic phase family
is preserved as a scoped theorem without a value gain; it has not been
promoted to a resolution of the catalyst question.

A possible, unpursued question suggested by the transition is whether
fixed-retention local moments for arbitrary symmetric Hadamard parents
admit a universal free-compression description. No such universality
or resulting cap curve was proved or used here. The sparse graph proof
alone must not be cited as a fixed-positive-retention moment limit.

No research-bearing files were created under temporary or ignored
paths during this closing audit. Existing canonical proofs, checkers,
and result JSON files remain in place for the root checkpoint.
