# Flatification campaign closing preservation / replay audit

2026-09-07. The 11:19:37 UTC snapshot preserves 1,891 research payloads,
155,772,363 bytes, all hashes verified. Later canonical proof/audit edits
are committed directly; no polishing was required for preservation.

The snapshot scan found no credential candidate or oversized unpreserved
research output. The same 14 historical/dynamic dependency flags remain
reviewed as in the preceding snapshots. The 176,489,665-byte raw witness
output remains locally available and byte-for-byte preserved by tracked
lossless gzip plus the reviewed original/compressed checksum manifest.

The final finite-certificate audit found an avoidable runtime dependency:
the four order-12 scripts named prebuilt binaries under tmp/. They now
build their tools from the tracked exact Gray-code and bridge evaluator
C++ sources via `flatify_construct_2026_09_07_build_tools.py`. Required
software is g++ with C++17, Python and NumPy from the documented venv.
Builds use a fresh directory under /home/math/quadra/tmp, never /tmp.
The new replay executable is excluded only with its pinned binary/source
hashes and explicit recipe in reviewed_build_products.json.

The finite phase certificate input is tracked JSON; the source signing is
tracked `extension_nested_m11_to_12.json`. No ignored mathematical input
is required for either exact verifier. The random search's elapsed-time
budget makes its output heuristic even with its recorded fixed RNG seed.
Its exact witnesses and exclusion certificates are independently replayable.
The director reran the complete phase verifier from a newly compiled
tracked source at closing: cap60 excluded for both polarities using
1092/1087 constraints; both cap62 witnesses passed all 8,388,608
projective assignments. The output exactly matches the committed JSON.

Research originals under tmp/ intentionally remain uncommitted in their
working locations; the tracked snapshots preserve them. Environments,
caches, reviewed software dependencies and regenerable products are
excluded. The original exact-range executable is also a reviewed build
product, not unpreserved mathematics. Archived source whitespace is left
unchanged deliberately; canonical new files passed the scoped diff check.
