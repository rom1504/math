# Publication audit: bounded polynomial-BH investigation

2026-09-16. Closing preservation audit. Mathematical conclusions are in
`artifacts/bh_2026_09_16_final_synthesis.md`, not inferred from this file.

## Preservation scope

The first checkpoint, commit `3bdaf9e`, was pushed normally. Its snapshot
`2026-09-16/060207Z` preserves1,877 payloads (154,816,421 bytes), with
all SHA-256 hashes verified. The closing snapshot also preserves later
proofs, agent drafts, primary literature, failed ideas and exact/diagnostic
outputs. Canonical proofs/scripts remain in their established locations.
The root's unproved closing escape checks are deliberately saved as a
scratch draft, not promoted to the theorem collection.

The pre-publication inventory found no credential candidate or oversized
research requiring a new storage decision. Environments/caches, reviewed
software dependencies and reviewed reproducible builds are excluded.
Preservation is not verification; old unfinished scripts are retained.

## New-code dependency review

- `bh_2026_09_16_power_checks.py`: NumPy and explicit tracked JSON
  witnesses at orders6--14. Exact caps are re-enumerated; optimality is
  not imported. All source paths and hashes occur in its saved output.
- `bh_2026_09_16_filter_checks.py`: NumPy, mpmath and tracked
  `computations/results/exact_m9.json`. Energy/Walsh/Chebyshev arithmetic
  uses Python integers; numerical coefficient norms are diagnostics.
- `bh_2026_09_16_orbit_flattening_checks.py` and
  `bh_2026_09_16_planted_tail_checks.py`: the tracked power-check module
  and its declared packages/witnesses; no ignored research input.
- `bh_mechanism_2026_09_16_joint_fields.py` and
  `bh_mechanism_2026_09_16_weighted_ports.py`: NumPy/SciPy;
  the weighted-port program additionally uses the tracked Walsh builder in
  `computations/principle_synthesis_2026_09_07_sparse_active_check.py`.
  Exact identities and floating diagnostics are separately labeled.
- `bh_mechanism_2026_09_16_chebyshev_recovery.py`: NumPy/SciPy and
  standard-library rational arithmetic. Its generic quadrature output
  is floating-point diagnostic data, not a quadratic-signing certificate.
- `bh_mechanism_2026_09_16_frame_saturation.py`: standard-library
  integer/rational arithmetic, deterministic seed and all selected
  rows/cores saved in the output. No external input.
- `bh_2026_09_16_independent_audit_checks.py`: standard-library exact
  arithmetic and the explicit independently re-enumerated H6 witness.

All nine replay programs pass Python compilation. The closing root replay
again passed36 random/Sylvester matrix checks, nine stored cap witnesses,
24 nonlinear filters, and the frame/orbit/planted/primary-contraction
checks. Exact certificates are not inferred from the77 floating norm
rows, five floating quadratures, or other numerical optimizations.
The primary paper PDFs/HTML and research notes are archived provenance,
not an execution dependency of the canonical replay programs.

The automatic ignored-reference scan is conservative and is not a proof
of dynamic dependency closure. New canonical imports and literal input
paths were manually reviewed as above. Historical unresolved references
remain documented in the existing archive audits and are not inputs to
the new mathematics.

## Intentionally retained working originals

- `computations/decisive_independent_exact_range_scan_2026_09_07`:
  reviewed compiled build, with source/hash provenance; not new research.
- `computations/results/flatify_director_rank_two_witness_profiles_2026_09_07_m32.json`:
  oversized raw research already preserved byte-for-byte in tracked
  lossless compression and its checksum manifest. The working original
  remains in place; it is not an unpreserved result.

No working original was deleted to make Git status appear clean.

## Closing verification

Snapshot `2026-09-16/070059Z` contains1,888 payloads and157,463,634 bytes;
all payload hashes verify. No concurrent-write warning was found. The71
ignored references comprise40 research,15 reviewed build references, two
directories/import paths and the same14 historical/dynamic flags. None
is a new missing input to the canonical BH replay programs.

The fresh, separately initialized refinement reviewer also completed a
standalone reconstruction and reran all31 exact checks without finding
a gap. Its artifact and the director's failed-escape scratch note are
included in the snapshot. Normal publication push and the post-commit
current-work audit are the remaining operational checks.
