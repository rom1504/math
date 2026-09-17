# Paper-combination campaign preservation audit

## First checkpoint

Snapshot `2026-09-17/172847Z` preserves1,911 research payloads, with
SHA-256 verification. Environments/caches, reviewed dependencies and
reproducible builds are excluded. The existing oversized raw experiment
remains covered by its byte-verified tracked compressed archive.

New build review: `tmp/paper_portfolio_2026_09_17/discrepancy/exact_fixed_signing_gray`
has SHA-256 `6611b17f46c5a59a85805c4e0fe73afb606762c762e527b9b8983a105643810a`;
tracked `computations/exact_fixed_signing_gray.cpp` has SHA-256
`5a08f4e89660efe23319351084d12007b4a44942fdcdc170d3c2fc05346c7586`.
The script rebuilds it with `g++ -O3 -std=c++17`. The unpublished snapshot
classification was corrected from unreviewed executable to build product;
only its redundant archived executable copy was removed, not the working
executable, source, inputs, logs or mathematical output.

New proof scripts have no ignored mathematical input. The isotropy program
reads tracked old witnesses and rebuilds the documented executable. Other
programs use explicit constants/enumeration and write their own JSON output.
Optional package versions are recorded in
`computations/paper_portfolio_requirements_2026_09_17.txt`.
The14 unresolved/dynamic dependency-reference flags are unchanged historical
entries, not new dependencies. Canonical proofs remain in artifacts/.

Remaining untracked working originals are retained intentionally: reviewed
compiled programs, the already losslessly archived oversized output, and
ongoing campaign work produced after the snapshot. They are not discarded
to make Git status empty. The campaign continues; the next substantive
checkpoint must preserve its later drafts and results again.
