#!/usr/bin/env python3
"""Audit one-edge cap-lowering repairs on exact exposed shells.

The repair criterion is (DR.22) in
``artifacts/cross_order_outward_director_review.md``.  Permutation
invariance lets us fix the restricted set to the first ``m`` vertices.
Switching quotient enumeration is exact through order eight.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from audit_canonical_disorder_root_gauge import (
    internal_edges,
    projected_masks,
    root_caps,
)


def internal_edge_toggles(n: int, m: int) -> list[int]:
    """Return root-gauge masks for flipping one edge inside ``[m]``."""
    positions = {edge: bit for bit, edge in enumerate(internal_edges(n))}
    output: list[int] = []
    for i in range(m):
        for j in range(i + 1, m):
            if i:
                output.append(1 << positions[(i, j)])
                continue
            # Flipping (0,j) and switching vertex j restores root gauge.
            toggle = 0
            for k in range(1, n):
                if k != j:
                    toggle ^= 1 << positions[tuple(sorted((j, k)))]
            output.append(toggle)
    return output


def audit(n: int, child_orders: list[int], levels: list[int]) -> dict:
    parent_caps = root_caps(n)
    records = []
    for m in child_orders:
        child_caps = root_caps(m)
        projection = projected_masks(n, m)
        toggles = internal_edge_toggles(n, m)
        for level in levels:
            shell = np.flatnonzero(parent_caps == level)
            normalized_level = level / n**1.5
            bad = shell[
                child_caps[projection[shell]] / m**1.5
                > normalized_level + 1e-12
            ]
            repairable = sum(
                any(parent_caps[int(mask) ^ toggle] <= level - 2
                    for toggle in toggles)
                for mask in bad
            )
            records.append(
                {
                    "N": n,
                    "m": m,
                    "cap_level": level,
                    "root_shell_size": int(len(shell)),
                    "bad_incidence_count": int(len(bad)),
                    "bad_incidence_mass": float(len(bad) / len(shell)),
                    "one_internal_edge_repairable_count": int(repairable),
                    "one_internal_edge_repairable_fraction": float(
                        repairable / max(1, len(bad))
                    ),
                }
            )
    return {
        "schema": "quadratic-signing-exposed-shell-edge-repair-audit-v1",
        "records": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--order", type=int, default=8)
    parser.add_argument("--children", type=int, nargs="+", default=[3, 4, 5])
    parser.add_argument("--levels", type=int, nargs="+", default=[10, 12])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rendered = json.dumps(
        audit(args.order, args.children, args.levels),
        indent=2,
        sort_keys=True,
    )
    if args.output:
        args.output.write_text(rendered + "\n", encoding="utf-8")
    else:
        print(rendered)


if __name__ == "__main__":
    main()
