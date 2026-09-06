# Fourth six-hour campaign: recursive construction and near-optimal seeds

Authorization: the user explicitly resumed six additional active hours from
the latest head b1745cc. Start2026-09-06 15:15:02 UTC; target end21:15:02 UTC.
The objective remains convergence or nonconvergence of M_n/n^(3/2), without
presuming a value. Checkpoints are steering events, not stopping points.

Writing this file successfully verifies workspace write access. Shell and
the repository Python environment are available; no execution blocker has
been found. Existing unrelated/untracked artifacts are preserved.

## Initial scope and verification

Read README.md, STEERING.md, the final synthesis and the strict all-order
upper theorem. Preserve the reported interval
[.4333221116640807,.499432220485404] and its proofs; independently reconstruct
the upper chain and correct any claim if that reconstruction exposes a gap.
Prior audit verdicts are not proof steps.

| Role | Initial target | Last result | Next falsifiable step |
|---|---|---|---|
| Reconstruction | Standalone proof of the new all-order upper theorem | Assigned | Rebuild dependencies and replay exact certificate |
| Seed transfer | All-order propagation of actual good finite signings | Assigned | Derive a seed-dependent construction and its precise cap loss |
| Adversarial/alternative | Test seed propagation and independent mechanisms | Assigned | Produce a quantitative actual-signing falsifier or a new implication |
| Director | Integrate the recursive factorization with original minimax data | Initial proof-obligation map | Identify where actual seed information enters or disappears |

The main unresolved implication is not a new numerical upper bound: a
construction must propagate seed caps approaching the original liminf to
all sufficiently large orders without a fixed normalized loss. No present
proof shows that the existing weave ensemble is asymptotically optimal.

## Checkpoint 1: reconstruction and exact limitations, resumed 17:02 UTC

The normal service hit a usage interruption around 16:05. The user said
continue; execution and all three researchers were available again at
17:02:38. This interval is not active research time. No model/account or
permission workaround was used. The approximate working end moves to 22:15.

The original interval is unchanged. The new upper theorem passed independent
reconstruction and the director's complete exact certificate replay. The
following are proved unless explicitly described as diagnostics:

| Result | Proof | Reproducible check |
|---|---|---|
| Standalone upper proof; direct-T stopping; H2/H12 order supply | [Reconstruction](transfer_reconstruction_standalone_2026_09_06.md) | `computations/transfer_reconstruction_exact_checks_2026_09_06.py --full-certificate` |
| Variable-precision representation and arbitrary orthogonal gate factorization | [Kernel/factorization](transfer_reconstruction_variance_kernel_factorization_2026_09_06.md), [Schur interface](transfer_seed_multivariate_schur_interface_2026_09_06.md) | `computations/transfer_seed_schur_twins_verify_2026_09_06.py` (Schur portion is floating diagnostic) |
| Arbitrary full sign base with variable row energy retained | [Tilted weave](transfer_seed_nonorthogonal_tilted_weave_2026_09_06.md) | Same seed checker, exact flatness portion |
| Near-optimal twins defeat all full-rank precision determinant bounds | [Proof](transfer_seed_twins_precision_obstruction_2026_09_06.md), [independent audit](transfer_reconstruction_twins_independent_audit_2026_09_06.md) | `computations/transfer_reconstruction_twins_exact_audit_2026_09_06.py` |
| Exact finite permanent and sqrt(15)/8 certificate floor | [Director proof](transfer_director_exact_permanent_floor_2026_09_06.md), [independent audit](transfer_reconstruction_permanent_floor_audit_2026_09_06.md) | `computations/transfer_reconstruction_permanent_floor_exact_checks_2026_09_06.py`; director checker separates rational from 80-digit diagnostics |
| Full/restricted weave and self-square actual seed losses | [Proof](transfer_adversary_seed_loss_2026_09_06.md) | `computations/transfer_adversary_seed_loss_2026_09_06.py` |
| Typical small restrictions of actual minimizers lose normalized cap | [Proof](transfer_adversary_random_restriction_2026_09_06.md) | `computations/transfer_adversary_random_restriction_2026_09_06.py` |

The director replayed every listed checker. Exact finite checks protect
arithmetic and conventions; the uniform and asymptotic claims have written
proofs. No new original lower or upper endpoint is claimed at this checkpoint.

| Track | Last concrete result | Next falsifiable step | Stop/redirect criterion |
|---|---|---|---|
| Reconstruction | Standalone theorem and independent floor/twin audits pass | Exact joint colored-edge contraction exponent | A lost leading term or infeasible repair invalidates the claimed sharpness |
| Seed transfer | Actual sign realization, but precision cost retains full Gram response | Random restrictions of literal seed tensors | A fixed lower witness gap or an unchanged full optimization obligation |
| Adversarial | R<=T mapped precisely; weighted 2-by-2 saturation imposes reciprocal spectra | Prove or falsify an asymptotic reciprocal array-pair construction | Finite low moments alone are not an asymptotic obstruction |
| Director | Exact underlying certificate floor, independent of recursive upper approximations | Uniform typical Gaussian spectrum theorem and actual annealed spin count | No inference from first moment to typical cap without another argument |

The exact original gap remains liminf seed landing with vanishing normalized
loss. A proposed next observable must enter a quantitative cap inequality.
This checkpoint is a steering event; substantive research continues.

## Preservation checkpoint: 2026-09-06 17:36 UTC

The user explicitly required durable preservation of research regardless of
temporary placement. README now requires archival preservation at substantive
checkpoints. The dated snapshot in
`research_archive/2026-09-06/173630Z/` preserves 1,886 research files with verified
SHA-256 hashes, including drafts, failed runs, seeds and subagent reports.
Canonical proofs/certificates were not moved. Two old unavailable chiral source
modules and other replay qualifications are documented in its dependency audit.
No endpoint claim is changed by this preservation operation. Agents continued
research during the inventory; later work will be preserved at the next
substantive checkpoint. This is not the end of the active campaign.

## Checkpoint 2: entropy-aware extraction and actual ensemble limits

The new rigorous all-order upper is <.499432211; the lower remains
.4333221116640807. This is an unconditional decimal improvement, NOT a proof
of convergence. The new ingredient is mandatory Gibbs-cluster entropy around
any quadratic extremizer; the original recursion and certificate are unchanged.

| Verified result | Proof | Reproducible checker |
|---|---|---|
| Gibbs-cluster upper extraction, exact signs/diagonal/depth/order limits | [Proof](transfer_reconstruction_gibbs_cluster_extraction_2026_09_06.md) | `transfer_reconstruction_gibbs_cluster_exact_checks_2026_09_06.py`; `transfer_adversary_cluster_cap_audit_2026_09_06.py` |
| Uniform typical Gaussian Hadamard spectra and actual annealed pressure | [Director](transfer_director_typical_hadamard_spectra_2026_09_06.md), [audit](transfer_reconstruction_typical_hadamard_pressure_audit_2026_09_06.md) | `transfer_reconstruction_typical_hadamard_exact_checks_2026_09_06.py` |
| Exact homogeneous endpoint pressure and hard-tail count | [Counting](transfer_reconstruction_homogeneous_edge_count_2026_09_06.md), [tail](transfer_reconstruction_signed_gaussian_tail_bridge_2026_09_06.md) | `transfer_reconstruction_homogeneous_edge_count_exact_checks_2026_09_06.py` |
| Near-identical Gaussian replicas falsify the proposed rate bound | [Proof](transfer_reconstruction_two_replica_rate_falsifier_2026_09_06.md) | `transfer_reconstruction_two_replica_rate_falsifier_2026_09_06.py` |
| Actual even self-tensors: complete classification, uniform weave-seed loss, positive-retention divergence | [Proof](transfer_seed_tensor_power_instability_2026_09_06.md) | `transfer_seed_schur_twins_verify_2026_09_06.py` |
| Extremely thin chosen tensors encode all signings; random ones forget seeds | [Proof](transfer_seed_very_thin_tensor_restrictions_2026_09_06.md) | `transfer_seed_thin_tensor_verify_2026_09_06.py` |

Checker paths are under `computations/`. The director independently replayed
all seven with `transfer_director_checkpoint2_replays_2026_09_06.py`; saved
output is `computations/results/transfer_director_checkpoint2_replays_2026_09_06.json`.
Exact and numerical portions remain labeled. The sharpness/count statements
are annealed, not lower bounds on typical signings. The old one-extremizer
floor does not apply to entropy-corrected extraction.

Next: local-field-adaptive entropy and conditional candidate exclusion,
intermediate tensor retention, and asymptotic nonlocal seed operators. These
newer tasks are pending director audit and are not smuggled into this table.
Convergence requires the same missing liminf-to-all-order comparison. Continue
through the active campaign, with new research archived at this checkpoint.

## Checkpoint 3: actual selector laws, deterministic escapes, correlation credit

The original interval remains [.4333221116640807,.499432211); convergence
and nonconvergence remain open. The following claims have written proofs,
independent reconstruction, and reproducible finite regressions where relevant.

| Verified result | Exact scope / proof |
|---|---|
| Every bounded-cap parent has sparse random restrictions of cap >=2/pi-o_P | [Actual-parent theorem](transfer_seed_sparse_bounded_cap_parent_theorem_2026_09_06.md), [director dependency audit](transfer_director_sparse_restriction_audit_2026_09_06.md) |
| Uniform small fixed-retention loss, probability <=K/n | [Connected-graph covariance proof](transfer_adversary_fixed_retention_random_loss_2026_09_06.md) |
| Good selectors have fraction <=3 exp(-kappa n^2/D); successful alternative laws pay the corresponding KL cost | [Convex-distance theorem](transfer_director_exponential_selector_cost_2026_09_06.md) |
| Any prescribed good sublinear child embeds in an asymptotically minimizing parent at O(n/D) normalized cost | [Finite planting lemma](transfer_adversary_deterministic_planting_escape_2026_09_06.md) |
| One Hadamard parent family has typical strict-subhalf 31/32 restrictions but sparse cap >=2/pi | [Retention transition](transfer_seed_random_retention_transition_2026_09_06.md) |
| Explicit even-spectral/local-diagonal defect forced by a low actual cap | [Finite covariance clipping](transfer_seed_exceptional_selector_finite_defects_2026_09_06.md) |
| Fixed-marginal clipped correlations give an extensive pressure gain; adaptive fields give additional product credit when L1 dispersion is present | [Correlation theorem](transfer_reconstruction_correlated_cluster_credit_2026_09_06.md), [adaptive fields and hub falsifier](transfer_reconstruction_adaptive_cluster_fields_2026_09_06.md) |
| Stationary centered-chaos identity; raw Gaussian logdet is not a lower bound | [Exact identity and scalable counterexamples](transfer_reconstruction_stationary_chaos_logdet_obstruction_2026_09_06.md) |
| Attenuated Gaussian-threshold pressure resummation with uniform O(sqrt N) error; unattenuated finite comparison false | [Threshold theorem](transfer_reconstruction_gaussian_threshold_resummation_2026_09_06.md) |
| Complement exposure can determine the whole weave cap; seed-only conditional tails need not have extensive speed | [Exposure](transfer_adversary_complement_exposure_2026_09_06.md), [seed-only falsifier](transfer_adversary_seed_only_conditional_tail_2026_09_06.md) |
| Two intermediate random self-tensor ranges have rigorous actual cap losses | [Range theorem](transfer_seed_intermediate_random_tensor_ranges_2026_09_06.md); no claim covering the intervening window |

All fourteen finite replays are preserved by
`computations/transfer_director_checkpoint3_replays_2026_09_06.py` in
`computations/results/transfer_director_checkpoint3_replays_2026_09_06.json`.
Numerical reconnaissance remains labeled numerical. The primary imported
graph and convex-distance theorems were checked at their exact hypotheses.
The director's unfinished free-projection and conditional-exposure calculations
are separately preserved as drafts, not silently converted into proof claims.

The meaningful distinction is now sharper: ordinary sampling cannot transfer
the optimizing coefficient, but specially correlated selections can be good,
even inside asymptotically minimizing parents. There is still no lossless
upward construction from a liminf seed. Correlated pressure credit is a positive
mechanism, but decimal gains alone are not that missing implication.

Continue the active six-hour campaign after this checkpoint. Next work tests
sharper discrete correlation entropy and an independently assessed original
convergence mechanism, rather than accumulating more selector diagnostics.

## Checkpoint 4: discrete information retained through quantization

The fixed-parameter discrete threshold-entropy target is now proved, with
independent reconstruction of its full operator argument and rational constants.
The order-uniform bound pays two rare-channel factors, not the larger continuous
Gaussian KL. It provides a strict further construction gain while leaving the
safely displayed interval [.4333221116640807,.499432211) unchanged.

| Verified result | Proof / audit |
|---|---|
| Threshold likelihood smoothing, positive operator comparison, random-restriction determinant bound, and useful coefficient 22 | [Theorem](transfer_seed_threshold_entropy_resampling_2026_09_06.md), [independent audit](transfer_adversary_threshold_entropy_resampling_audit_2026_09_06.md) |
| Finite quantizers and arbitrary squared-entry geometry satisfy a weighted information inequality | [Director theorem](transfer_director_quantized_gaussian_information_2026_09_06.md), [audit](transfer_adversary_quantized_information_audit_2026_09_06.md) |
| Exact integration into the old annealed upper certificate, with subextensive errors and unchanged limit orders | [Integration](transfer_director_threshold_cap_integration_2026_09_06.md) |
| Every fixed retention has some liminf-realizing near-minimizers with positive-probability leading child loss | [First-exit theorem](transfer_seed_all_fixed_retention_nearmin_loss_2026_09_06.md), [audit](transfer_adversary_all_fixed_retention_audit_2026_09_06.md); not an exact-minimizer or every-seed claim |
| All Paley conferences have isotropic absolute-ground laws; optimized one-row extensions cost at least sqrt(n), parity-rounded | [Fresh proof](transfer_fresh_cavity_and_isotropic_ground_law_2026_09_06.md); no asymptotic optimality assumption |

The director replayed five programs with
`computations/transfer_director_checkpoint4_replays_2026_09_06.py`.
The saved result separates rational assertions, exact integer enumeration,
and numerical diagnostics. No replay has an ignored research input.

The rare-threshold paired example only excludes a multiplier uniform as
delta tends to zero; the earlier broader interpretation was corrected. A
full Hadamard family satisfies the new information theorem but still has
cap tending to 1/2, demonstrating that a separate annealed upper budget is
essential. This theorem does not silently become a universal cap upper bound.

The preservation audit additionally identified one source-less executable,
`tmp/search_r33`; it is retained as research. The other 26 executables have
reviewed source mappings and explicit build recipes. The archiver now pins
binary and source hashes instead of excluding every ELF automatically.

Continue the active campaign. The fresh lower cross-order attempt and positive
seed-transfer work remain separate from the completed entropy subproblem.

## Checkpoint 5: exact minimizing quantifiers, balanced selection, and local traps

The rigorous asymptotic interval is unchanged and convergence remains open.
The following results have independently reconstructed proofs; finite screens
and imported solver statuses are kept separate from those proofs.

| Result | Scope and proof |
|---|---|
| Orientation-balanced near-minimizers at EVERY order, Q<=M_n+O(n^(11/8)), abs(P-R)=O(n^(11/8)) | [Director construction](transfer_director_balanced_orientation_repair_2026_09_06.md), [independent audit](transfer_seed_orientation_repair_independent_audit_2026_09_06.md); no isotropy or op-norm claim |
| Actual near-minimizers can forbid isotropy on every o(n^(4/3))-slack shell, with mean-slack and covariance lower bounds | [Clique theorem](transfer_adversary_nearmin_clique_orientation_gap_2026_09_06.md), [independent audit](transfer_seed_nearmin_clique_independent_audit_2026_09_06.md) |
| Powered pointwise aggregation fails at linear scale for strict-subhalf parent and children, the small child EXACTLY minimizing | [Stopped asymmetric theorem](transfer_fresh_subhalf_asymmetric_power_obstruction_2026_09_06.md), [audit](transfer_adversary_stopped_asymmetric_power_audit_2026_09_06.md); not an M-value or comparable-split counterexample |
| Balanced pointwise aggregation still fails after truncating child caps at a fixed subhalf coefficient | [Capped theorem](transfer_fresh_balanced_capped_power_obstruction_2026_09_06.md); no claim that the actual child caps are subhalf |
| Orientation gap n is compatible with strict single-edge and sub-sqrt(n) Hamming local optimality | [Hadamard trap](transfer_adversary_exact_minimizer_orientation_gap_2026_09_06.md); these are NOT global minimizers |
| Exact-ground isotropy fails for every n4 minimizer; exact rational finite response LPs at fifteen stored cases | [Elementary proof and certificate scope](transfer_adversary_exact_minimizer_isotropy_finite_2026_09_06.md) |

The original-value reverse-Fekete criterion and its all-order proof are in
`transfer_fresh_reverse_fekete_2026_09_06.md`. The exact-parent lexicographic
partition attempt proves genuine exact-ground blockers but no simultaneous
large-block control; see `transfer_fresh_exact_minimum_partition_attempt_2026_09_06.md`.
The finite [partition screen](transfer_seed_existential_partition_screen_2026_09_06.md)
passes its powered error tests but does not prove a growing-order assertion.
The stronger demand for both children to be exact already fails for one stored
order-12 minimizer at the comparable split four/eight.

The director's full-constraint orientation-gap MILP reported infeasibility
at n9 and n10, and timed out at n11--14 without a witness. Those statuses are
exploratory only, not independent infeasibility certificates. The code,
parameters and informative outputs are retained regardless of that failure.

An [attributed external audit](transfer_adversary_external_finite_frontier_audit_2026_09_06.md)
independently replays exact Q15=27 and Q16=30 witnesses. It does not import
M15=27 or M16>=28: full catalogue completeness remains external. Pinned URLs,
source hashes and that trust boundary are recorded. Our mathematical data and
verifiers are preserved; unlicensed downloaded dependencies remain locally with
per-file pinned source manifests rather than being republished.

The director replayed all seven checkpoint programs with
`computations/transfer_director_checkpoint5_replays_2026_09_06.py`;
the complete arguments and outputs are saved in its matching results JSON.
The external witness verifier was additionally rerun in its standalone,
preserved-data mode, checking both full projective histograms and every
single-vertex deletion without reading the external cache.

An independently checked continuation result, recorded separately in
`transfer_adversary_isotropic_response_extension_2026_09_06.md`, proves
`0<=I(B)-I(A)<=r sqrt(n)+(r/2)sqrt(r-1)` for ANY r-vertex extension,
where `I` is maximum absolute energy averaged under an isotropic Boolean law.
This prevents sublinear vertex extensions from erasing a hypothetical
macroscopic isotropic-response gap. It does not assert such a gap for
minimizers. Its finite checker has 2046 comparisons and 76 exact rational
primal/dual values.

Continue the active campaign after this checkpoint. Actual-feedback covariance
and the precisely defined catalyst comparison are the next bounded positive
and adversarial tests. The improved description of failed partition arguments
is not counted as a proved original convergence step.

## Checkpoint 6: same-order balance and a genuinely signed return

The original interval remains unchanged, and convergence is OPEN. Two
positive structural obligations were removed without being called convergence.

| Verified result | Exact scope and dependencies |
|---|---|
| Same-order Q<=M_n+O(n^(5/4)) and orientation gap O(n) | [Director repair](transfer_director_same_order_orientation_repair_2026_09_06.md), [independent finite/iterative audit](transfer_adversary_same_order_orientation_repair_audit_2026_09_06.md) |
| Balanced actual near-minimizers can still require isotropic mean slack Omega(n^(4/3)) | [Scoped counterexample](transfer_adversary_balanced_nearmin_isotropic_defect_2026_09_06.md); not exact minimizers and not a leading n^(3/2) gap |
| Inhomogeneous Gaussian responses retain a nonlinear return despite mixed row-coefficient signs | [Trace-rank theorem](transfer_director_inhomogeneous_gaussian_return_rigidity_2026_09_06.md), [independent audit](transfer_adversary_gaussian_return_rigidity_audit_2026_09_06.md) |
| Actual first marked feedback has liminf[Q(A)/(n sqrt(n-1))-j_n]>=gamma(f,L)>0 | [Actual mixed and joint comparison](transfer_seed_high_degree_actual_feedback_probe_2026_09_06.md), [independent reconstruction](transfer_adversary_high_degree_actual_feedback_audit_2026_09_06.md), [endpoint implication](transfer_fresh_conditional_actual_feedback_endpoint_gain_2026_09_06.md) |

For the last theorem f is a fixed odd ternary function of the original
first marked Gaussian fields, with nonzero occupied and hole masses; the
actual signing has fixed normalized operator cap L. A fixed high-original-
degree channel and finite odd-Hermite catalog give a positive FULL BC
correlation. A signed tensor of rank at most n controls the aggregate before
absolute values. The proof then checks joint independence of the probe from
the literal regressed query list, rather than inferring it from marginal
Gaussianity. Both endpoint orientations have a positive instability density;
exact thinned flips give the cap gain. Finite tie-induced means are retained.

The director read and reconstructed both the mixed-contraction audit and
the independently derived quantitative endpoint proof. Fixed source, degree,
variance cutoff and root caps precede the order limit. This theorem does NOT
yet cover the richer 200-degree lower-certificate frame, and its conservative
gain does not establish domination of the O(1/L) principal-core loss.

Six programs are replayed by
`computations/transfer_director_checkpoint6_replays_2026_09_06.py`, with
complete output in the matching results JSON. Exact integer/rational tests
are separated from floating Gaussian and rank diagnostics. The director's
new frozen f(Y) test exhausts 33,272 seed configurations on fifteen actual
stored signings, checking the two endpoint orientations and the ascent
identity exactly; twelve have positive finite ascent in both orientations.
This is a finite regression, not proof of the asymptotic comparison.

The catalyst pass did not settle R=T. Its exact odd-Walsh saturation test
and finite inverse-spectrum enumeration are preserved, with primary-source
normalizations and clear finite-versus-asymptotic boundaries. The exact-
minimizer isotropy dual attempt is likewise preserved as unsuccessful.

Continue after this checkpoint, approximately through 22:15 UTC. The next
discriminating tests concern richer-frame bounded comparison and compact-
source uniformity, not tuning the old terminal functional or restarting
the failed finite catalyst search. Research-bearing drafts and failures
are archived again while this work continues.
