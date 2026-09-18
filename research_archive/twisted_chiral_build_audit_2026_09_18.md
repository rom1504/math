# Twisted chiral campaign: build and dependency audit

2026-09-18. The build-product candidate report in
computations/results/twisted_chiral_2026_09_18_build_products_candidate.json
records compiler version, exact source and executable hashes, and the full
build command for the campaign's C++ programs. Sources are tracked. The
director reviewed the cap evaluator, root-gauge census formula and drivers.
Search drivers rebuild missing/stale executables; scratch executables are
not undocumented mathematical input. Frozen source copies retain the exact
algorithms used by already-running computations. Frozen-run diagnostics are
research provenance and are preserved separately unless byte-identical
generated copies have been reviewed.

The six reviewed binaries are search, search_frozen,
search_bound_initial_frozen, verify, classify9 and classify10 under
tmp/twisted_chiral_2026_09_18/. Their source/binary hashes are registered in
reviewed_build_products.json. This is a reproducibility review, not a proof
of the mathematical conclusions of every program.

All canonical Python calculations read explicit constants or tracked
matrices. Source filenames and seed indices are retained in computational
outputs. The two director inverse JSONL files originate in scratch and must
be included in the dated preservation archive. They are generated outputs,
not ignored prerequisites of the inverse-audit source program.

The pre-existing untracked September7 executable remains a reviewed build
product. The pre-existing 176489665-byte raw profile JSON remains covered by
its tracked byte-verified lossless archive and manifest. No originals are
deleted to clean the worktree.
