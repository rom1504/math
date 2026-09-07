# Fourth active-campaign preservation checkpoint

1,905 research payloads, 176,506,678 bytes. This is preservation, not an
upgrade of any draft's evidentiary status. Canonical proofs remain in place.
All payload hashes are verified with preserve.py --verify. Subsequent
canonical additions are committed directly; later drafts enter the next
snapshot. The same 14 historical/dynamic dependency flags remain; the new
standalone certificates use tracked scripts and the recorded Python packages.

One additional oversized original, the order-2048 heuristic witness JSON,
is preserved losslessly in tracked gzip form. The manifest category
preserved_compressed_research identifies it; reviewed_compressed_research.json
pins both hashes. The preservation utility checks exact decompressed bytes.
It is NOT excluded as a disposable build product. The original stays local.

Environments, caches, downloaded software dependencies and source/hash-reviewed
executables are excluded. No credential-pattern candidate was detected.
