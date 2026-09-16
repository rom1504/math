# BH-paper investigation: first preservation checkpoint

UTC2026-09-16 06:02:07. This is preservation, not a blanket verification
of every archived claim. The two-hour mathematical investigation remains
active until07:13:45UTC.

All1,877 research payload hashes were independently replayed with
`research_archive/preserve.py --snapshot research_archive/2026-09-16/060207Z --verify`.
The snapshot includes154,816,421 bytes of proofs, ongoing drafts, failed
experiments, downloaded primary papers, exploratory code and outputs.
Canonical BH proofs and replay scripts are committed in their existing
artifact/computation locations; later edits will be preserved again.

No credential candidate or unpreserved oversized research output was
found. The old176,489,665-byte raw rank-two output remains byte-verified
in its tracked lossless gzip and manifest. The old exact-range executable
remains an excluded build product with its tracked source/build audit.
Neither working original was deleted.

The conservative dependency scan has the same14 historical/dynamic
references previously reviewed. New BH scripts use standard-library,
numpy/scipy/mpmath, tracked exact witness JSON, and the documented tracked
Walsh builder; they do not require an ignored research input. Primary
downloads and independent replay outputs placed under tmp/ are retained
here as provenance, not silently required at replay time.

Inventory, original paths, SHA-256 hashes, exclusions, and environment
versions are in the adjacent manifest/runtime files. Do not infer proved
status from archival placement; see each canonical proof's audit labels.
