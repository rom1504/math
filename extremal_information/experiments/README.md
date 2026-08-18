# Exact finite laboratory

All programs in this directory are deterministic and use the repository-local
Python environment.  They are falsification tools; no finite pattern is
silently extrapolated.

Run from the repository root:

```bash
.venv/bin/python extremal_information/experiments/entropy_overlap_lab.py
.venv/bin/python extremal_information/experiments/pinned_query_rate_verify.py
.venv/bin/python extremal_information/experiments/build_quadratic_landscape_dataset.py
.venv/bin/python extremal_information/experiments/verify_code_replica_hierarchy.py
.venv/bin/python extremal_information/experiments/verify_phase2_response_geometry.py
.venv/bin/python extremal_information/experiments/verify_phase2_code_syndrome_profiles.py
.venv/bin/python extremal_information/experiments/verify_phase2_normalized_code_rate_distortion.py
.venv/bin/python extremal_information/experiments/verify_phase3_closed_algebra_claims.py
.venv/bin/python extremal_information/experiments/verify_phase3_prime_cycle_summary.py
.venv/bin/python extremal_information/experiments/verify_phase3_tropical_defect_saturation.py
.venv/bin/python extremal_information/experiments/verify_phase3_geodesic_fibre_bound.py
.venv/bin/python extremal_information/experiments/verify_phase3_geodesic_second_audit.py
.venv/bin/python extremal_information/experiments/verify_phase3_matroid_quotients.py
.venv/bin/python extremal_information/experiments/verify_invariant_rational_grid_shadowing.py
.venv/bin/python extremal_information/experiments/verify_exposed_commuting_germ_lower_bound.py
.venv/bin/python extremal_information/experiments/verify_reward_congruence_nonlattice.py
.venv/bin/python extremal_information/experiments/verify_pwa_binary_counter_exposure.py
.venv/bin/python extremal_information/experiments/verify_reward_congruence_coloring.py
.venv/bin/python extremal_information/experiments/verify_word_spectrum_pathlift_gap.py
```

The scripts write only to this directory by default.

## `entropy_overlap_lab.py`

This verifies three independent claims exactly.

1. The codes
   `C={0000,0011,0101,0110}` and
   `D={0000,0011,0101,1001}` have the same ordered pair-distance enumerator
   and covering radii two and three.  Their Cartesian-power enumerators and
   additive radius separation are checked symbolically.
2. Two saved order-eight complete signings have identical energy histograms
   but caps `16` and `20` after the same all-negative one-vertex extension.
3. Every residual unlabeled graph in the NetworkX graph atlas is checked for
   orders four through eight.  No equal exact
   `(energy(x),energy(y),overlap(x,y))` signature has two different multisets
   of one-vertex extension caps in this range.

Output:
[`entropy_overlap_results.json`](entropy_overlap_results.json).

The third statement is a complete finite census only through order eight.

## `pinned_query_rate_verify.py`

This enumerates all `2^binom(n,2)` sign interactions through order five.  For
every landscape and every field direction it verifies the pinning identity,
recovers each edge coefficient by its degree-two Walsh coefficient, and
checks that all response vectors are distinct.

Output:
[`pinned_query_rate_results.json`](pinned_query_rate_results.json).

This script is a finite falsifier for the algebra in the information theorem;
the rate lower bound itself is proved analytically.

## `build_quadratic_landscape_dataset.py`

This builds the order-eight experimental dataset.  Switching makes the first
row positive; the remaining negative edges are represented by each unlabeled
graph on seven vertices from the NetworkX graph atlas.  The script groups all
`1044` residual graph representatives into `243` exact
energy--energy--overlap signature classes.

For one representative of every class it records:

- the exact matrix and energy histogram;
- positive, negative, and absolute extrema;
- exact near-cap counts;
- trace powers through degree eight;
- row sums; and
- the full histogram of one-vertex extension caps.

Output:
[`quadratic_landscape_order8.json`](quadratic_landscape_order8.json).

The dataset is intended for collision search and theorem discovery.  It is
not a machine-learning certificate and the chosen representative is not a
canonical switching/permutation normal form.

## `verify_code_replica_hierarchy.py`

This supplies exact finite checks for Theorem 3.3.  It verifies every proper
selected-codeword column profile for the parity constructions at
`r=3,5,7,9`, exhaustively computes covering radii at `r=3,5`, checks the
alternating-binomial identity, and enumerates the full ambient `t`-point
membership/distance census for the four-bit base pair.  The latter agrees
through `t=4` and first differs at `t=5`.

Output:
[`code_replica_hierarchy_results.json`](code_replica_hierarchy_results.json).

## `verify_phase2_response_geometry.py`

This exhaustively checks the finite constants and counterexamples used in the
second theory-building phase:

- inverse-Hamming response moduli for boundary kernels, arbitrary rooted
  codes, shifted quadratics, and Max-Cut;
- the equal-outer-spectrum code collision;
- the rare-matching conditional-variance calculation;
- same-space response cancellation; and
- the `D_r` tropical-rank versus mean-square scale separation.

Output:
[`phase2_response_geometry_results.json`](phase2_response_geometry_results.json).

## `verify_phase2_code_syndrome_profiles.py`

This enumerates binary syndrome profiles for small interface dimensions,
checks min-plus convolution against support union, finds two length-five
nonisometric fragments with the same syndrome profile, and verifies the
equal-outer-spectrum pair separated by the full-rank environment used in
Example 6.

Output:
[`phase2_code_syndrome_profiles_results.json`](phase2_code_syndrome_profiles_results.json).

## `verify_phase2_normalized_code_rate_distortion.py`

This exhaustively checks the block-direct-sum construction in Theorem 8.3 for
five small `(L,q)` pairs.  It verifies full rank and common padded length, the
exact subset-count covering-radius formula, the induced directed-difference
response metric, and its inverse-Hamming lower bound.  In total it checks 104
state--environment pairs and 42 unordered state pairs.

Output:
[`phase2_normalized_code_rate_distortion_results.json`](phase2_normalized_code_rate_distortion_results.json).

## `verify_phase3_closed_algebra_claims.py`

This exhaustively checks two phase-three algebra claims on small binary
instances.  First, it compares the proposed closed formula for the complete
future residual-rank response metric with every pair of flats and every
future flat through width four.  Second, it checks the carrier-span
compression pointwise and after maximizing over targets, for one to three
carriers, every background dictionary, and every future dictionary through
width three.

Output:
[`phase3_closed_algebra_claims_results.json`](phase3_closed_algebra_claims_results.json).

## `verify_phase3_tropical_defect_saturation.py`

This checks finite-valued min-plus kernel powers on small cyclic groups,
distinguishes raw Hamming-ball smoothing error from the additional
one-blur-to-many-blur defect, and verifies the fixed-chart syndrome quotient
on all support pairs through width three and a deterministic width-four
sample.  It also records why the identically-infinite profile must be excluded
from the finite sup-norm statement.

Output:
[`phase3_tropical_defect_saturation_results.json`](phase3_tropical_defect_saturation_results.json).

## `verify_phase3_prime_cycle_summary.py`

This checks the finite phase meshes and the exact closed-algebra obstruction
for prime cyclic groups.  It records the actual uniform response-net error and
the optimal one-state closed error for several primes and distortions.

Output:
[`phase3_prime_cycle_summary_results.json`](phase3_prime_cycle_summary_results.json).

## `verify_phase3_geodesic_fibre_bound.py`

This exhaustively checks every spanning binary support through width four and
every diametral shortest representation.  It verifies independence, the exact
zero fibre, diameter-two nonzero fibres, and the `D+1` anticode bound for
`D>=3`; it also records the sharp `D=2` exception and deterministic random
checks at widths five and six.

Output:
[`phase3_geodesic_fibre_bound_results.json`](phase3_geodesic_fibre_bound_results.json).

## `verify_phase3_geodesic_second_audit.py`

This independent implementation enumerates small Hamming anticodes, checks
the hard-core quotient against every raw (including nonspanning) future at
width three, and compares the actual width-four state count with the proved
enumeration bound.

Output:
[`phase3_geodesic_second_audit_results.json`](phase3_geodesic_second_audit_results.json).

## `verify_phase3_matroid_quotients.py`

This verifies fixed-subspace and triggered contraction fibres and exhaustively
enumerates all `3,616` join congruences of the binary width-three subspace
lattice.  It checks the canonical kernel factorization, pullback recovery, and
the exact residual-rank oscillation formula class by class.

Output:
[`phase3_matroid_quotients_results.json`](phase3_matroid_quotients_results.json).

## `verify_phase3_geodesic_synchronization.py`

Checks the binary geodesic/cycle criterion through width three, the sharp
fibre-stripping family through quotient diameter four, all-future response
comparisons on sampled contexts, and the coupled-bent construction at
`(D,k)=(6,4)`.

Output:
[`phase3_geodesic_synchronization_results.json`](phase3_geodesic_synchronization_results.json).

## `verify_phase3_vector_blr_response_audit.py`

Independently checks scalar and vector BLR constants, all cycle-contracting
maps on `F_2^3` modulo linear coordinates, dense and partial all-future word
profiles through ambient width three, and the sharp `2h` selector family.

Output:
[`phase3_vector_blr_response_audit_results.json`](phase3_vector_blr_response_audit_results.json).

## `verify_phase3_geodesic_sync_full_audit.py`

This independent audit verifies the exact selector-cube rooted-response
formula through quotient dimension four.  It also checks the coupled-bent
geodesic construction at `k=4` and at the first block extension `k=8`,
including all additive triangles, exact distance to linear sections, and the
full Cayley diameter.

Output:
[`phase3_geodesic_sync_full_audit_results.json`](phase3_geodesic_sync_full_audit_results.json).

## `verify_phase3_transversal_composition.py`

Checks the constant-error synchronization replacement, affine-state and
affine-rank profile bounds, the exact mixed-cycle identity, the equal-depth
radius obstruction through `D=12`, and the two linear-code pair-spectrum
examples.

Output:
[`phase3_transversal_composition_results.json`](phase3_transversal_composition_results.json).

## `verify_phase3_transversal_composition_audit.py`

Independently enumerates affine feature-state products and joins, tests their
all-future midpoint decoders, and checks the equal-state macroscopic
obstruction and exact support-separating queries.

Output:
[`phase3_transversal_composition_audit_results.json`](phase3_transversal_composition_audit_results.json).

## `verify_phase3_mixed_circuit_hierarchy.py`

Checks antipodal elimination, the circuit and nullity bounds, mixed-holonomy
rank under random fragment partitions, the strict arity hierarchy through
seven sources, sharp disjoint-circuit accumulation, and a finite one-channel
rooted-response packing.

Output:
[`phase3_mixed_circuit_hierarchy_results.json`](phase3_mixed_circuit_hierarchy_results.json).

## `verify_phase3_multichannel_holonomy_packing.py`

Exhaustively checks the parallel-pair response formula and its
Hamming--Hausdorff comparison on all two-dimensional subspaces of
`F_2^5`.  It then verifies the macroscopic witness separation for all `155`
two-subspaces of the binary Reed--Muller `[16,5,8]` host code.

Output:
[`phase3_multichannel_holonomy_packing_results.json`](phase3_multichannel_holonomy_packing_results.json).

## `verify_phase3_carrier_capacity.py`

Exhaustively verifies the carrier/Hausdorff response comparison for all `242`
zero--one weighted carriers of the five-cycle.  It also constructs the
field-multiplication subspace of `4 x 4` binary matrices and checks the
rank-metric carrier packing on all `35` of its two-dimensional subspaces.

Output:
[`phase3_carrier_capacity_results.json`](phase3_carrier_capacity_results.json).

## `verify_phase3_carrier_relation_response_law.py`

Checks the metric-carrier rough isometry and its two sharp boundaries on
finite-field examples: exact collapse for a redundant discrete alphabet, a
two-state quotient for a linear-diameter two-scale metric, and positive Lee
and flag-ultrametric host packings.

Output:
[`phase3_carrier_relation_response_law_results.json`](phase3_carrier_relation_response_law_results.json).

## `verify_phase3_metric_quotient_synchronization.py`

Exhausts all `728` zero--one presented carriers on a six-point product
metric, all `65,535` nonempty carriers in the `2 x 2` binary rank-metric row
projection, the two-scale finite-field quotient, and `104,976` min-plus
nonamplification instances.

Output:
[`phase3_metric_quotient_synchronization_results.json`](phase3_metric_quotient_synchronization_results.json).

## `verify_phase3_scale_rank_response_sandwich.py`

Enumerates every subspace of the three-dimensional binary two-scale carrier
to recover its complete scale-rank function.  It also constructs the binary
`[4 x 4,8,3]` Gabidulin rank-metric host, checks all `255` one-dimensional
carrier profiles and their `32,385` pairs, and attains the predicted response
gap.

Output:
[`phase3_scale_rank_response_sandwich_results.json`](phase3_scale_rank_response_sandwich_results.json).

## `verify_phase3_scale_rank_sandwich.py`

Independently checks the strict-threshold Singleton obstruction, the complete
two-scale rank curve and projected decoder, small Hamming separated ranks,
and separate sharp examples for the fibre, lift, and presentation terms.

Output:
[`phase3_scale_rank_sandwich_results.json`](phase3_scale_rank_sandwich_results.json).

## `verify_phase3_qary_multichannel_holonomy.py`

Exhaustively checks over `F_3` that scalar-closed quotient fibres introduce no
extra word-length shortcuts: all `624` independent ordered pairs in
`F_3^3`, all `16,848` kernel-endpoint profiles, and the response--Hausdorff
comparison for all `78` pairs of projective lines are verified.

Output:
[`phase3_qary_multichannel_holonomy_results.json`](phase3_qary_multichannel_holonomy_results.json).

## `verify_phase3_multichannel_response_entropy.py`

Independently checks the exact weighted-Cayley formula, the classification of
exact profile collisions by distinct heavy generator columns, failure of
general channel-basis invariance, the span sandwich, and `74,924` direct
two-fragment endpoint identities on small binary carriers.

Output:
[`phase3_multichannel_response_entropy_results.json`](phase3_multichannel_response_entropy_results.json).

## `verify_phase3_sparse_flat_grassmannian.py`

Checks the exact sparse-flat directed-ball identity at every center and
integer threshold for `(D,k)=(4,2),(5,2)`.  It also verifies two centers with
isometric quotient leader geometry but different symmetric balls, and the
scalable rooted-kernel diameter gap.

Output:
[`phase3_sparse_flat_grassmannian.json`](phase3_sparse_flat_grassmannian.json).

## `verify_phase3_hamming_grassmannian_coding_barrier.py`

Exhaustively verifies the binary line metric, exact line-ball formula, and
`A_2-1 <= P <= A_2` through `D=5`, plus `4,152` systematic-chart comparisons.

## `verify_phase3_hamming_grassmannian_injection_barrier.py`

Checks `25,650` subspace pairs for the low-weight sum-code certificate,
the finite common-host counterexample, monotonicity, and `62,807` numerical
regressions of the analytic entropy domination theorem.

Output:
[`phase3_hamming_grassmannian_injection_barrier_results.json`](phase3_hamming_grassmannian_injection_barrier_results.json).

## `verify_phase3_hamming_grassmannian_falsifier.py`

Verifies the line slice through `D=6`, all directed comparisons in the
seven-plane, simplex-line, and rank-multiplication alphabets, and length-two
product additivity.

Output:
[`phase3_hamming_grassmannian_falsifier_results.json`](phase3_hamming_grassmannian_falsifier_results.json).

## `verify_phase3_finite_alphabet_response_amplification.py`

Independently checks the exact directed response tables for seven Hamming
simplex lines and seven `F_8` multiplication lines, together with all `1,176`
length-two word pairs in each model.

Output:
[`phase3_finite_alphabet_response_amplification.json`](phase3_finite_alphabet_response_amplification.json).

## Dynamic selector checks

`verify_invariant_rational_grid_shadowing.py` checks 120,066 exact
grid/word inequalities for two switched rational clamps. It verifies exact
forward invariance and no growth of the initial shadow error despite
slope-one cells.

`verify_exposed_commuting_germ_lower_bound.py` checks both the block-product
prototype and the stronger fixed-three-letter construction. Two permutation
selectors generate every constant-weight binary state, and a repeatable
centered-coordinate probe separates every pair at reward rate one.

`verify_reward_congruence_nonlattice.py` exhaustively checks the two
incomparable bounded-error partitions of the reset automaton through depth
12, the bounded pairwise raw-response distance, and the linear drift of their
one-block join.

`verify_pwa_binary_counter_exposure.py` checks the dual-rail increment map
through ten bits, verifies that every cyclic probe shift is distinct, and
tests projective homogeneity and span-one invariance on small rational grids.
The theorem proves the construction at every width.

`verify_reward_congruence_coloring.py` exhausts all 1,099 graphs through five
vertices and 54,253 graph/partition pairs. It verifies that the identity-
dynamics reward construction is feasible at error `1/2` exactly on
independent-set partitions, so its minimum quotient size equals chromatic
number.

`verify_word_spectrum_pathlift_gap.py` checks 15,024 finite word products for
the exact scalar-spectrum/full-row-lift separation and 6,216 periodic-window
and unique-row claims for the fixed-binary de Bruijn family. It verifies that
wordwise critical witnesses need not assemble into one reusable path lift.
It also verifies that every tested length-`2m` product has one all-zero row
and every other row identically `-m`, hence projective rank one despite the
full rooted path-lift lower bound.

`verify_zero_relation_survival_carrier.py` exhausts every two-letter
two-state zero-relation alphabet and checks equivalence between nonmortality
of endpoint subsets and cyclicity of every relation product. It also verifies
the fixed-two-letter `2^r` mortality-monitor construction through `r=8` and
the one-state support carrier for binary de Bruijn systems through memory
seven.

`verify_ising_anticipatory_support.py` checks the strict width-two Ising
benchmark: two backward-surjective supports, no one-support realization,
four distinct forward row signatures, and all 1,092 words through length six
for both the flat response and the interacting formula
`2N_ca^cyc-N_a`.

## `actual_child_orbit_posterior_quotient.py`

Exhaustively evaluates all `2^16` bridges for each of the two certified
order-eight pressure-minimizer classes paired with the order-two child at
raw temperature `t=3` and negative exponent `lambda=1`. Integer group and
profile enumeration gives `19` and `22` simultaneous signed-similarity
orbits. A modular certificate over `F_1000003`, followed by transcendence of
`e^3`, proves that the exact averaged posterior distinguishes every orbit.
The same computation proves that the `12` combined-energy shells are too
coarse. This is an exact finite result, not an all-order orbit bound.

Output:
[`../../computations/results/actual_child_orbit_posterior_quotient.json`](../../computations/results/actual_child_orbit_posterior_quotient.json).

## `actual_child_orbit_scaling.py`

Extends the exact posterior-orbit audit to every certified thermal pressure
minimizer available from orders two through eight at raw `t=3`.  It combines
an exact cap-gap/energy-histogram pressure certificate with signed-group,
rooted-profile, and modular posterior enumeration.  The posterior saturates
the simultaneous symmetry quotient in every class except order three,
where an exact denominator-type cancellation proves one all-temperature
collision.  Order nine is recorded only as an explicitly ineligible cap-
minimizer diagnostic.

Output:
[`../../computations/results/actual_child_orbit_scaling.json`](../../computations/results/actual_child_orbit_scaling.json).

## `actual_child_escort_low_degree_falsifier.py`

Exhaustively selects the thermal-pressure-minimizing children at physical
raw temperature `4/sqrt(N)` for `3+3`, `4+4`, and `3+7`, enumerates both
orientations of the complete bridge cube, and solves the weighted normal
equations for the best degree-one and degree-three polynomial approximation
to every edge-orbit cavity under the `lambda=1` inverse escort.  It also
checks the exact parity law: even Walsh degrees are orthogonal to every
cavity response.  The residuals are complete-cube numerical evidence, not
an interval certificate or asymptotic lower bound.

Run:

```bash
.venv/bin/python \
  extremal_information/experiments/actual_child_escort_low_degree_falsifier.py
```

Output:
[`../../computations/results/actual_child_escort_low_degree_falsifier.json`](../../computations/results/actual_child_escort_low_degree_falsifier.json).

## `actual_child_posterior_frame_wind_tunnel.py`

Reconstructs the exact edge posterior-mean matrix from deleted-edge cavity
responses for the certified `beta=4` thermal-minimizing children at splits
`3+3`, `4+4`, and `3+7`, with `lambda=1` negative-escort weighting over the
complete bridge cube.  It tests fixed child spectral frames, elementary
`B1`/`B^T1` field frames, and the exact natural signed-symmetry quotient.
Individual posterior matrices are often strongly rank one, but all three
proposed frame synchronizations fail at the finite frontier.  Exact Burnside
orbit counts are `18`, `1438`, and `16148`; SVD-derived quantities are
numerical finite evidence, not asymptotic claims.

Run:

```bash
.venv/bin/python \
  extremal_information/experiments/actual_child_posterior_frame_wind_tunnel.py
```

Audit:
[`actual_child_posterior_frame_wind_tunnel_report.md`](actual_child_posterior_frame_wind_tunnel_report.md).

Output:
[`../../computations/results/actual_child_posterior_frame_wind_tunnel.json`](../../computations/results/actual_child_posterior_frame_wind_tunnel.json).
