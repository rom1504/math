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
