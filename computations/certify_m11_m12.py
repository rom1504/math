#!/usr/bin/env python3
"""Assemble and verify the campaign certificates for M_11 and M_12.

The lower bound for M_11 is solver-certified: the symmetry-complete cap-15
CP-SAT model is infeasible.  The upper bound is an explicit cap-17 matrix.
For M_12, monotonicity gives M_12 >= M_11, all order-12 energies are even,
and an explicit one-vertex extension has cap 18.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from exact_mn_milp import exact_profile, stable_matrix_hash


def matrix_from(payload: dict[str, object]) -> np.ndarray:
    for key in ("matrix", "parent_matrix"):
        if key in payload:
            return np.asarray(payload[key], dtype=np.int8)
    raise KeyError("matrix or parent_matrix")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--m11-lower",
        type=Path,
        default=Path("computations/results/cpsat_m11_cap15.json"),
    )
    parser.add_argument(
        "--m11-upper",
        type=Path,
        default=Path("computations/results/nested_10_in_11_cap17.json"),
    )
    parser.add_argument(
        "--m12-upper",
        type=Path,
        default=Path("computations/results/extension_nested_m11_to_12.json"),
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    lower = json.loads(args.m11_lower.read_text())
    upper_11 = json.loads(args.m11_upper.read_text())
    upper_12 = json.loads(args.m12_upper.read_text())
    if not (
        lower.get("n") == 11
        and lower.get("decision_cap") == 15
        and lower.get("solver", {}).get("status") == "INFEASIBLE"
        and lower.get("model", {}).get("root_gauge") is True
        and lower.get("model", {}).get("basic_permutation_and_complement_symmetry")
        is True
    ):
        raise AssertionError("unexpected M_11 lower-bound certificate")

    matrix_11 = matrix_from(upper_11)
    matrix_12 = matrix_from(upper_12)
    profile_11 = exact_profile(matrix_11)
    profile_12 = exact_profile(matrix_12)
    if len(matrix_11) != 11 or profile_11["M"] != 17:
        raise AssertionError((len(matrix_11), profile_11["M"]))
    if len(matrix_12) != 12 or profile_12["M"] != 18:
        raise AssertionError((len(matrix_12), profile_12["M"]))
    if not np.array_equal(matrix_12[:11, :11], matrix_11):
        raise AssertionError("the M_12 witness is not the claimed M_11 extension")
    if stable_matrix_hash(matrix_11) != upper_11["matrix_sha256"]:
        raise AssertionError("M_11 hash mismatch")
    if stable_matrix_hash(matrix_12) != upper_12["parent_matrix_sha256"]:
        raise AssertionError("M_12 hash mismatch")

    payload = {
        "schema": "quadratic-signing-certified-values-m11-m12-v1",
        "classification": "solver-certified exact values with exhaustively verified upper witnesses; no standalone lower-bound proof object",
        "values": {"11": 17, "12": 18},
        "M11_implication": [
            "the symmetry-complete cap-15 CP-SAT model is INFEASIBLE",
            "the explicit order-11 witness has exhaustively recomputed cap 17",
            "therefore M_11=17",
        ],
        "M12_implication": [
            "deleting a vertex and maximizing over its spin proves M_(n+1)>=M_n",
            "therefore M_12>=M_11=17",
            "all order-12 energies have the parity of binom(12,2)=66 and hence are even",
            "the explicit order-12 extension has exhaustively recomputed cap 18",
            "therefore M_12=18",
        ],
        "sources": {
            "M11_lower": str(args.m11_lower),
            "M11_upper": str(args.m11_upper),
            "M12_upper": str(args.m12_upper),
        },
        "witnesses": {
            "11": {
                "matrix_sha256": stable_matrix_hash(matrix_11),
                "profile": profile_11,
            },
            "12": {
                "matrix_sha256": stable_matrix_hash(matrix_12),
                "profile": profile_12,
            },
        },
    }
    print("verified solver-certified exact values M_11=17 and M_12=18")
    print(f"M_11 witness {payload['witnesses']['11']['matrix_sha256']}")
    print(f"M_12 witness {payload['witnesses']['12']['matrix_sha256']}")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
