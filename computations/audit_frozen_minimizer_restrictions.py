#!/usr/bin/env python3
"""Audit restrictions of the exact uniform minimizing fibre at orders 7--8.

The stored orbit classification contains every root-gauged labelled
minimizer.  Reintroducing all switching gauges and restricting to the first
``m`` vertices gives the exact marginal of a uniformly selected labelled
minimizer.  This script measures how much mass leaves the smaller minimizing
fibre and the frozen terms in the canonical-disorder restriction identity.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from analyze_equal_split_partitions import add_orbit_class  # noqa: E402
from audit_canonical_disorder_restriction import caps, edges, entropy  # noqa: E402


def root_minimizers(payload: dict) -> list[np.ndarray]:
    n = int(payload["order"])
    classes: list[dict[str, set[bytes]]] = []
    for row in payload["classes"]:
        add_orbit_class(classes, np.asarray(row["representative_matrix"], dtype=np.int8))
    keys: set[bytes] = set()
    for orbit in classes:
        keys.update(orbit["plus"])
        keys.update(orbit["minus"])
    expected = int(payload["minimizing_signing_count"])
    if len(keys) != expected:
        raise AssertionError((len(keys), expected))
    return [np.frombuffer(key, dtype=np.int8).reshape(n, n) for key in keys]


def restriction_mask(matrix: np.ndarray) -> int:
    mask = 0
    for bit, (i, j) in enumerate(edges(len(matrix))):
        if matrix[i, j] < 0:
            mask |= 1 << bit
    return mask


def audit(source: Path, max_child: int) -> list[dict[str, float | int]]:
    payload = json.loads(source.read_text(encoding="utf-8"))
    n = int(payload["order"])
    parent_cap = int(payload["target_cap"])
    parents = root_minimizers(payload)
    rows: list[dict[str, float | int]] = []
    first_m = max(2, math.ceil(n / 3))
    last_m = min(n - 1, max_child, math.floor(2 * n / 3))
    for m in range(first_m, last_m + 1):
        child_caps = caps(m)
        child_min = int(np.min(child_caps))
        counts = np.zeros(len(child_caps), dtype=np.int64)
        for parent in parents:
            restricted = parent[:m, :m]
            for switch_mask in range(1 << (m - 1)):
                switch = np.ones(m, dtype=np.int8)
                for i in range(1, m):
                    if (switch_mask >> (i - 1)) & 1:
                        switch[i] = -1
                switched = switch[:, None] * restricted * switch[None, :]
                counts[restriction_mask(switched)] += 1
        marginal = counts / float(np.sum(counts))
        expected_cap = float(np.dot(marginal, child_caps))
        nonminimal_mass = float(np.sum(marginal[child_caps > child_min]))
        q = (m * (m - 1)) / (n * (n - 1))
        parent_entropy = math.log(len(parents) * (1 << (n - 1)))
        shearer_slack = entropy(marginal) / q - parent_entropy
        energy_excess = (
            math.sqrt(m) * expected_cap / q - math.sqrt(n) * parent_cap
        )
        rows.append(
            {
                "N": n,
                "m": m,
                "parent_minimizer_count_root_gauge": len(parents),
                "parent_cap": parent_cap,
                "child_cap": child_min,
                "nonminimal_restriction_mass": nonminimal_mass,
                "expected_restriction_cap": expected_cap,
                "restriction_entropy": entropy(marginal),
                "shearer_slack": shearer_slack,
                "frozen_energy_excess": energy_excess,
            }
        )
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "sources",
        nargs="*",
        type=Path,
        default=[
            HERE / "results" / "m7_minimizer_orbits.json",
            HERE / "results" / "m8_minimizer_orbits.json",
        ],
    )
    parser.add_argument("--max-child", type=int, default=6)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rows = [row for source in args.sources for row in audit(source, args.max_child)]
    rendered = json.dumps({"records": rows}, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(rendered + "\n", encoding="utf-8")
    else:
        print(rendered)


if __name__ == "__main__":
    main()
