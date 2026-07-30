#!/usr/bin/env python3
"""Certify the conference completion of the saved order-13 cap-20 bridge.

If a signing A of order 13 satisfies 13I-A^2=ss^T for a sign vector s and
As=0, then

    C = [[0, s^T], [s, A]]

is a symmetric conference matrix of order 14.  This script verifies those
identities in exact integer arithmetic and exhaustively evaluates the Boolean
caps of A, C, and every one-vertex deletion of C.
"""

from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path

import numpy as np

from exact_mn_milp import exact_profile, stable_matrix_hash


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = json.loads(args.source.read_text())
    key = "parent_matrix" if "parent_matrix" in payload else "matrix"
    a = np.asarray(payload[key], dtype=np.int64)
    n = len(a)
    if n != 13:
        raise ValueError("this certificate expects an order-13 signing")
    defect = n * np.eye(n, dtype=np.int64) - a @ a
    if np.any(np.diag(defect) != 1):
        raise AssertionError("spectral defect does not have unit diagonal")
    s = defect[:, 0].copy()
    if not np.array_equal(defect, np.outer(s, s)):
        raise AssertionError("spectral defect is not a sign outer product")
    if not np.array_equal(a @ s, np.zeros(n, dtype=np.int64)):
        raise AssertionError("kernel-vector identity failed")
    conference = np.zeros((n + 1, n + 1), dtype=np.int64)
    conference[0, 1:] = conference[1:, 0] = s
    conference[1:, 1:] = a
    if not np.array_equal(
        conference @ conference, n * np.eye(n + 1, dtype=np.int64)
    ):
        raise AssertionError("conference square identity failed")
    profile_a = exact_profile(a)
    profile_c = exact_profile(conference)
    deletion_caps = []
    for vertex in range(n + 1):
        keep = [i for i in range(n + 1) if i != vertex]
        deletion_caps.append(exact_profile(conference[np.ix_(keep, keep)])["M"])
    exact_partitions = []
    for subset in combinations(range(n + 1), 7):
        if 0 not in subset:
            continue
        complement = tuple(i for i in range(n + 1) if i not in subset)
        left = conference[np.ix_(subset, subset)]
        right = conference[np.ix_(complement, complement)]
        if exact_profile(left)["M"] == 9 and exact_profile(right)["M"] == 9:
            exact_partitions.append((subset, complement))
    if not exact_partitions:
        raise AssertionError("no exact 7+7 child partition found")
    subset, complement = exact_partitions[0]
    left = conference[np.ix_(subset, subset)]
    right = conference[np.ix_(complement, complement)]
    bridge = conference[np.ix_(subset, complement)]
    output = {
        "schema": "quadratic-signing-conference-completion-m13-v1",
        "classification": (
            "proved exact integer conference completion and exhaustive finite "
            "Boolean profiles"
        ),
        "source": str(args.source),
        "order_13_matrix_sha256": stable_matrix_hash(a),
        "kernel_sign_vector": [int(value) for value in s],
        "identities": [
            "13I-A^2=ss^T",
            "As=0",
            "C=[[0,s^T],[s,A]] and C^2=13I",
        ],
        "order_13_profile": profile_a,
        "conference_matrix": [
            [int(value) for value in row] for row in conference
        ],
        "conference_matrix_sha256": stable_matrix_hash(conference),
        "order_14_profile": profile_c,
        "order_13_deletion_caps": deletion_caps,
        "exact_7_7_unordered_partition_count": len(exact_partitions),
        "first_exact_7_7_partition": {
            "left_vertices": list(subset),
            "right_vertices": list(complement),
            "left_matrix": [[int(value) for value in row] for row in left],
            "right_matrix": [[int(value) for value in row] for row in right],
            "bridge": [[int(value) for value in row] for row in bridge],
            "child_caps": [9, 9],
            "parent_cap": 21,
            "two_thirds_defect": 21 ** (2.0 / 3.0) - 2 * 9 ** (2.0 / 3.0),
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"A13 cap={profile_a['M']} C14 cap={profile_c['M']} "
        f"deletion_caps={deletion_caps} exact_7+7={len(exact_partitions)}"
    )
    print(
        f"A13 hash={output['order_13_matrix_sha256']} "
        f"C14 hash={output['conference_matrix_sha256']}"
    )
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
