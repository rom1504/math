#!/usr/bin/env python3
"""Exhaustively audit the universal block doubling on saved child signings.

For every symmetric zero-diagonal sign matrix S, not only conference matrices,

    D(S) = [[S, S+I], [S+I, -S]]

is a valid signing of twice the order.  Its exact energy can be written as

    2 H_S(x) - 4 H_(S[J])(x_J) + n - 2|J|.

This script computes the exact doubled caps for supplied children when the
parent is small enough for exhaustive enumeration.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from conference_double_construction import double_conference
from exact_mn_milp import exact_profile, stable_matrix_hash


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--parent-order-limit", type=int, default=20)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    records = []
    for path in args.inputs:
        payload = json.loads(path.read_text())
        matrix = np.asarray(payload["matrix"], dtype=np.int8)
        child = exact_profile(matrix)
        parent = double_conference(matrix)
        if len(parent) > args.parent_order_limit:
            continue
        parent_profile = exact_profile(parent)
        defect = parent_profile["M"] ** (2.0 / 3.0) - 2 * child["M"] ** (
            2.0 / 3.0
        )
        record = {
            "source": str(path),
            "child_order": len(matrix),
            "child_M": child["M"],
            "parent_order": len(parent),
            "parent_M": parent_profile["M"],
            "parent_normalized_M": parent_profile["M"] / len(parent) ** 1.5,
            "two_thirds_defect": defect,
            "parent_matrix_sha256": stable_matrix_hash(parent),
        }
        records.append(record)
        print(
            f"{record['child_order']}->{record['parent_order']}: "
            f"M {record['child_M']}->{record['parent_M']} "
            f"normalized={record['parent_normalized_M']:.12f} "
            f"defect={defect:+.12f}"
        )
    output = {
        "schema": "quadratic-signing-universal-double-audit-v1",
        "classification": "proved exhaustive finite caps for the saved child representatives",
        "exact_energy_identity": (
            "H_D(x,J)=2 H_S(x)-4 H_(S[J])(x_J)+n-2|J|"
        ),
        "records": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
