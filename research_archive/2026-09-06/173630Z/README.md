# Preservation checkpoint — 2026-09-06 17:36:30 UTC

**Status: archival snapshots, not newly verified mathematics.** Parent research
commit: `018675a`. Files were copied non-destructively while the three research
agents continued. No concurrent-copy warnings occurred. All 1,886 saved
payloads passed size and SHA-256 verification after copying.

## Contents and exclusions

- **1,886 research files / 153,016,506 bytes:** scratch derivations, exploratory
  Python/C++ programs, unsuccessful searches, finite outputs, seeds/witnesses,
  old and current subagent reports, literature provenance, and working proof
  drafts. Original paths are under `files/`; use the [manifest](manifest.json)
  to locate a file without guessing its status.
- Large groups include the August retrieval panel and its reports,
  actual-child experiments, exact-minimizer tail reruns, row certificates,
  Grothendieck/AMP source audits, and September recursive-weave/seed-transfer
  explorations. Failed and superseded runs were retained, not deduplicated by
  mathematical interpretation.
- **17,594 environment/cache files / 577,858,019 bytes excluded.** This includes
  `.venv`, the old PDF-tool environment, pip cache, and Python caches.
- **77 build products / 1,463,672 bytes excluded:** compiled executables and
  other reproducible compiled/package products. Their inventory remains in the
  manifest. Eight environment symlinks were recorded but not followed.
- No strong credential-pattern candidates and no research file above 25 MiB
  were found. Thus no oversized research was deferred to external storage.
  This scan is not a guarantee that arbitrary secrets can be detected.

The historical source files are not upgraded to current proof dependencies by
being archived. Some active proofs will change after this snapshot; subsequent
substantive checkpoints must preserve those changes too.

## Dependency audit

The automated [reference inventory](ignored_references.json) found 42 literal
`tmp/` references in committed Markdown/code: 23 refer to preserved research,
12 to compiled build products, two to preserved directories, and five needed
manual classification. All five are resolved as follows:

- `tmp/rich_core_optimized_replay.json` is an output path in a documented
  replay command, not a required input. The policy and canonical certificate
  are tracked.
- `tmp/tree_pair_audit` (two references) and `tmp/power16` are example compiled
  executables, reproducible from the C++ source containing the commands.
- `tmp/rp_audit.json` was explicitly removed in the old radial-path audit.
  The audit says it reproduced the tracked
  `computations/results/actual_child_radial_path_curvature_falsifier.json`.
  Its historical bytes cannot be recovered from this workspace; the canonical
  result and experiment source remain tracked. We have not invented a replay
  and represented it as the missing historical file.

An independent read-only Python/import and C++/input review additionally found:

| Finding | Resolution / remaining limitation |
|---|---|
| Two committed probes import untracked `continued_convergence_universal_overlap_2026_09_06.py` and `continued_convergence_latent_supersolution_test_2026_09_06.py` | Both sources preserved here and added at their existing `computations/` paths in this checkpoint. Their experimental status is unchanged. |
| `artifacts/package_chiral_equality_instances.py`, `search_chiral_equality_constraints.py`, and `verify_chiral_equality_bounded_checkpoint.py` require `audit_dependent_4lift` and/or `bounded_chiral_equality_search` | Both modules are absent from the workspace, including ignored files; no matching path appears in `git log --all`. These old scripts remain blocked for replay. Their NPZ witnesses are tracked. This is a documented historical preservation gap, not a new mathematical disproof. |
| `continued_convergence_restricted_weave_2026_09_06.py` runs a `tmp/` executable | It rebuilds the binary from tracked `continued_convergence_weave_ascent_2026_09_06.cpp`. Other old evaluator drivers similarly accept a compiled executable built from their cited tracked C++ sources. |
| `fresh_hadamard36_symmetrization_probe.py` downloads external matrix text | Its source URL is explicit in that script; replay depends on that external source remaining available. No corresponding local source cache was found. This old experiment is not an input to either current endpoint. |
| `requirements.txt` lists only NumPy/SciPy/OR-Tools, while some older experiments also use mpmath, cvxpy, SymPy, and NetworkX | [runtime.json](runtime.json) records all installed package names/versions without credentials or environment variables. It is preservation metadata, not a claim that every experiment has been reproduced. |

### Current endpoint dependency closure

The latest strict upper proof chain and its exact certificate have no essential
ignored-file input. The exact certificate imports the tracked ternary-latent
interval module and NumPy; its canonical result is tracked. Its extra `tmp/`
JSON files are independent replay outputs, now also preserved.

The reported `0.4333221116640807` lower certificate has an eight-script local
dependency closure, all tracked. Its optimized rectangle policy is tracked at
`computations/results/resumed_response_rich_core_optimized_rectangle_policy_2026_09_06.json`.
The exact replay target for that reported decimal is
`4333221116640807/10000000000000000` with `--theta 1`; the older `4333/10000`
example checks only a weaker value. The canonical optimized-policy certificate
is tracked. The analytic proof-dependency files for both endpoints were also
checked for presence/tracking. This preservation audit is **not** a fresh proof
of either endpoint; mathematical reconstruction belongs in canonical audits.

## What remains outside this commit

Working originals in `tmp/` and other untracked locations are deliberately left
in place. Their checkpoint versions are durable in this archive; active files
may subsequently differ. Canonical in-progress proofs/scripts are not promoted
to verified mathematics merely to clean status. Environments, caches,
credentials (if any), and rebuildable binaries remain excluded. The two missing
chiral sources cannot be committed because they are unavailable.

No guarantee is made for undiscovered dynamic dependencies: the exact scope of
the automated and independent audits is stated above. Future checkpoints must
record any additional dependency discovered during replay.
