# Research preservation archive

This is an append-only preservation layer, not a second proof ledger. A file
being saved here does **not** make its claims verified. Snapshots deliberately
include drafts, abandoned ideas, failed computations, raw reports, numerical
outputs, parameters, seeds, and literature used during research. Canonical
proofs, certificates, and scripts remain in their original locations.

## Index

| UTC checkpoint | Scope | Status / audit |
|---|---|---|
| [2026-09-07 06:05:13](2026-09-07/060513Z/README.md) | 1,868 research files; 154,968,179 bytes; renewed flatification work, improved ternary certificate, independent replays, MUB drafts and failed experiments | All payload hashes verified; new dynamic filename flag manually resolved; no credential candidate or oversized research output. Campaign continues. |
| [2026-09-07 04:32:20](2026-09-07/043220Z/README.md) | 1,847 research files; 153,117,802 bytes; final ramp/ETF work, affine-Gaussian gap, exact checks, drafts, failures and outputs | Payloads verified; all remaining research matched the current-work audit. New proof dependencies manually checked; no credential candidate or oversized output. |
| [2026-09-07 04:02:20](2026-09-07/040220Z/README.md) | 1,839 research files; 153,090,851 bytes; bipartite lower theorem, independent upper replay, midpoint ramp surgery, exact tests and continuing drafts | Preservation only; payload hashes verified. No credential candidate or oversized research output. Later edits will be preserved at the next checkpoint. |
| [2026-09-07 03:36:00](2026-09-07/033600Z/README.md) | 1,853 research files; 153,150,768 bytes; exact log-cosh identities, fresh lower replay, entropy surgery, universality falsifier, sparse benchmark and ongoing drafts | Preservation only; all payload hashes verified. No credential candidate or oversized research output. Canonical proof/certificate paths remain unchanged. |
| [2026-09-07 02:37:17](2026-09-07/023717Z/README.md) | 1,848 research files; 152,919,876 bytes; quadratic fluctuation proofs, width-pressure comparison, soft-flatness reduction, independent audits and continuing failed calculations | Preservation only; canonical status remains explicit. New compiled scan excluded by pinned source/binary hash; no oversized research or credential candidate. |
| [2026-09-07 01:44:16](2026-09-07/014416Z/README.md) | 1,859 research files; 152,954,617 bytes; exact optimizer counterexamples, interpolation proofs, unsuccessful polarity solver logs and continuing drafts | Preservation only; canonical proof status remains explicit. No oversized research or credential candidate. |
| [2026-09-07 00:43:21](2026-09-07/004321Z/README.md) | 1,864 research files; 153,049,957 bytes; exact Gaussian phase, transport, counterexamples, code, raw outputs and ongoing drafts | Preservation only; canonical status remains in proofs. Environments, caches and reviewed build/dependency files excluded. |
| [2026-09-06 23:39:40](2026-09-06/233940Z/README.md) | 1,932 research files; 155,077,507 bytes; new precision proof, exact certificates, independent reconstruction, failed seed comparisons and ongoing drafts | Preservation only; canonical mathematical status is explicit in each artifact. Environments and reviewed dependencies excluded by existing manifest rules. |
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

Latest publication audit: [2026-09-07 closing report](publication_audit_2026_09_07.md).
All 1,824 remaining uncommitted research originals matched saved payload hashes
after the closing mathematical checkpoints were pushed. Working originals
remain in place intentionally; the tracked archive is their durable copy.

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
