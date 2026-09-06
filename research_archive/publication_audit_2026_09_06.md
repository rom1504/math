# Closing preservation publication audit

2026-09-06, approximately 22:15 UTC. The substantive checkpoint was committed
and pushed as `c2bb2b92ed84fab67e386dca212870b3807d9a3d`. Local HEAD and
`origin/main` agreed. This report is published in the subsequent closing commit.

Seven dated snapshots preserve successive research states. The final snapshot,
`2026-09-06/220601Z`, contains 1,915 research payloads, 154,814,524 bytes,
with all payload hashes verified. Two later subagent audit amendments were
committed canonically rather than overwriting historical snapshot payloads.
The archive does not confer verified status on unfinished or failed work.

Post-push `preserve.py --audit-current` reported:

```json
{
  "matching_current_research": 1880,
  "not_in_snapshot": [],
  "changed_since_snapshot": [],
  "credential_candidates_needing_review": []
}
```

## What remains uncommitted and why

- 1,880 research working originals, including ignored `tmp/` files and
  untracked artifacts/outputs: retained for ongoing work, byte-identical to
  tracked archival payloads. They were not deleted or moved merely to clean
  Git status. Canonical proofs and certificates remain in their existing paths.
- 17,597 environment/cache files: excluded as nonresearch runtime material.
- 76 reviewed compiled/bytecode products: excluded with source/build
  provenance. The executable `tmp/search_r33`, whose source is unavailable,
  is preserved as research rather than treated as reproducible.
- 25 reviewed downloaded dependency files: locally retained, excluded from
  republishing, with exact hashes and pinned-commit retrieval locations in
  `reviewed_external_dependencies.json`. Research notes, independent verifiers
  and informative outputs are not covered by that exclusion.
- Eight symlinks: all belong to `.venv/` or the old PDF-extraction environment;
  their targets are recorded. No research symlink target was omitted.

No oversized research output required separate durable storage in this
inventory. No credential candidate was detected; the pattern scan is not
a guarantee against every possible secret format.

## Dependency and evidence boundaries

The literal ignored-reference audit and manual current-proof/replay review
found no undocumented ignored research input in the current theorem chains.
The six final replay programs use documented Python packages and tracked
source inputs. Their recorded source hashes and evidence labels distinguish
exact integer/rational checks from floating numerical diagnostics.

Two historical modules, `audit_dependent_4lift` and
`bounded_chiral_equality_search`, were already absent from the workspace and
available Git history. Their dependent old scripts remain blocked for replay;
this is documented in the first snapshot, not silently repaired or declared
reproducible. Neither module enters the current rigorous interval. External
catalogue completeness and arbitrary dynamic imports were not certified.

The README now requires this preservation sweep at substantive checkpoints
and before cleanup. The mathematical campaign continued alongside these
snapshots; its conclusions and remaining convergence gap are in the tracked
[final synthesis](../artifacts/transfer_campaign_final_synthesis_2026_09_06.md).
