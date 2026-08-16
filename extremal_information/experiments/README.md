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
