# Research preservation archive

This is an append-only preservation layer, not a second proof ledger. A file
being saved here does **not** make its claims verified. Snapshots deliberately
include drafts, abandoned ideas, failed computations, raw reports, numerical
outputs, parameters, seeds, and literature used during research. Canonical
proofs, certificates, and scripts remain in their original locations.

## Index

| UTC checkpoint | Scope | Status / audit |
|---|---|---|
| [2026-09-06 17:36:30](2026-09-06/173630Z/README.md) | 1,886 research files; 153,016,506 bytes; all available untracked/ignored research plus modified tracked files | Preservation only; SHA-256 verified. Two unavailable old source modules documented, not reconstructed. |
| [2026-09-06 17:55:33](2026-09-06/175533Z/README.md) | 1,898 research files; checkpoint-2 proofs/replays and newer agent drafts included | Preservation only; SHA-256 verified. Canonical verification status is in the campaign index, not inferred from this snapshot. |
| [2026-09-06 18:51:23](2026-09-06/185123Z/README.md) | 1,915 research files; 154,762,575 bytes; checkpoint-3 proofs, failed comparisons, numerical reconnaissance and unfinished calculations | Preservation only; all hashes verified. New canonical replay programs have no ignored research input. |
| [2026-09-06 19:29:59](2026-09-06/192959Z/README.md) | 1,901 research files; 154,677,428 bytes; threshold-information proofs, scoped counterexamples, failed reconnaissance, and the source-less experiment binary | Preservation only; all hashes verified. Executable exclusions now require reviewed binary and source hashes. |
| [2026-09-06 20:35:22](2026-09-06/203522Z/README.md) | 1,916 research files; 155,235,485 bytes; all-order orientation repair, exact-minimizer tests, failed MILPs, subagent drafts and independent finite replays | Preservation only; all hashes verified. External downloaded dependencies have per-file pinned-commit manifests; independent order-15/16 witness replay is standalone. |
| [2026-09-06 21:14:00](2026-09-06/211400Z/README.md) | 1,903 research files; 154,918,023 bytes; same-order repair, actual signed-feedback proofs, exact finite checks, abandoned catalyst tests and continuing drafts | Preservation only; all hashes verified. New mathematical and diagnostic scripts have tracked inputs; no new external dependency. |
| [2026-09-06 22:06:01](2026-09-06/220601Z/README.md) | 1,915 research files; 154,814,524 bytes; cap-only energy approximation, rich-feedback analytic lemmas, exact falsifiers, all pending drafts and failed searches | Preservation only; payload hashes and remaining current-work hashes checked. Six standalone replays have documented package/ tracked-source inputs. |

Each snapshot has `manifest.json` (original paths, hashes, byte counts,
provenance category, exclusions), `files/` (unchanged payloads),
`ignored_references.json` (conservative dependency-reference scan), and
`runtime.json` (observed package versions, not a tested lockfile).

Research literature/source downloads are retained as research provenance;
software environments, downloaded software dependencies, caches, credentials,
and regenerable compiled executables are excluded. No research payload in the
initial inventory exceeds the 25 MiB per-file review threshold. Oversized
research must have a checksum manifest and independently durable storage
location before it can be omitted from a snapshot. An ignored local path is
not such a location.

The [build-product audit](build_product_audit_2026_09_06.md) maps all 27
existing research executables and 50 bytecode files to preserved or vendor
sources. One executable, `tmp/search_r33`, has no recovered source and is
preserved as research. ELF format alone is no longer grounds for exclusion:
`reviewed_build_products.json` pins both the binary and reviewed source hashes.
New or changed executables are retained until that provenance review is done.
Build recipes are documented separately from claims of historical bit-for-bit
reproducibility; no working executable was deleted or overwritten.

Reviewed, byte-identical external software/source dependencies can instead be
preserved by a pinned upstream commit, durable retrieval URL and SHA-256 in
`reviewed_external_dependencies.json`. This applies to the inspected external
order-15/16 certificate package, whose upstream license is unspecified; it does
not exclude our derivations, independently authored verifiers, witness data or
informative results. Matching is per-file and per-hash, never a blanket exclusion
of a research directory. Local downloads are retained, not deleted. The external
certificate audit records the remaining completeness and replay boundaries.

## Checkpoint procedure

Run from the repository root with the project venv:

```bash
.venv/bin/python research_archive/preserve.py
.venv/bin/python research_archive/preserve.py --save
.venv/bin/python research_archive/preserve.py --snapshot research_archive/DATE/TIMEZ --verify
.venv/bin/python research_archive/preserve.py --snapshot research_archive/DATE/TIMEZ --audit-current
```

Review the inventory and exclusions before publishing. The utility refuses to
save detected credential candidates or oversized research needing a durable
storage decision; its pattern checks are not a complete secret audit. Review
concurrent-write warnings and resnapshot changed work. Add an index entry and
brief dependency findings, then commit and push. Do not wait for a draft to be
polished, and do not remove working originals to clean Git status.
The current-work audit reports new or changed uncommitted research not captured
by the selected snapshot. Run it after publication as well: already committed
canonical files need no duplicate working-copy preservation, while remaining
untracked/ignored research must still match a saved payload.

To restore one missing ignored dependency without overwriting newer work:

```bash
.venv/bin/python research_archive/preserve.py \
  --snapshot research_archive/2026-09-06/173630Z \
  --restore-path tmp/transfer_reconstruction_exact_certificate_2026_09_06.json
```

The `--restore-path` value is the original path in the manifest. Relative
imports and data paths in archived scripts were not rewritten; restore their
recorded dependencies when replaying them. Some archival experiments are
unfinished or broken by design. Consult the snapshot audit before relying on
them. A literal-path scan cannot establish full dependency closure for dynamic
code; important theorem chains require a manual dependency audit as well.
