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

## Checkpoint procedure

Run from the repository root with the project venv:

```bash
.venv/bin/python research_archive/preserve.py
.venv/bin/python research_archive/preserve.py --save
.venv/bin/python research_archive/preserve.py --snapshot research_archive/DATE/TIMEZ --verify
```

Review the inventory and exclusions before publishing. The utility refuses to
save detected credential candidates or oversized research needing a durable
storage decision; its pattern checks are not a complete secret audit. Review
concurrent-write warnings and resnapshot changed work. Add an index entry and
brief dependency findings, then commit and push. Do not wait for a draft to be
polished, and do not remove working originals to clean Git status.

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
