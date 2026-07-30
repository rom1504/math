#!/usr/bin/env python3
"""Assemble the symmetry-split CP-SAT certificate for M_13 and M_14.

After rooting an order-13 signing, global negation followed by re-gauging
complements its internal graph on 12 vertices.  One may therefore choose a
representative whose minimum internal negative degree lies in {0,...,5}.
The six supplied CP-SAT decisions must certify cap 18 infeasible in each
degree case.  Together with the explicit cap-20 bridge this proves M_13=20;
monotonicity, odd order-14 energy parity, and the explicit cap-21 conference
completion then prove M_14=21.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from exact_mn_milp import exact_profile, stable_matrix_hash


def matrix_from_payload(payload: dict[str, object]) -> np.ndarray:
    for key in ("matrix", "parent_matrix", "conference_matrix"):
        if key in payload:
            return np.asarray(payload[key], dtype=np.int8)
    raise KeyError("no signing matrix in payload")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("decisions", nargs=6, type=Path)
    parser.add_argument("--m13-witness", type=Path, required=True)
    parser.add_argument("--m14-witness", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    cases = []
    seen = set()
    for path in args.decisions:
        payload = json.loads(path.read_text())
        degree = payload["model"]["root_negative_degree_case"]
        if payload["n"] != 13 or payload["decision_cap"] != 18:
            raise AssertionError(f"wrong decision instance in {path}")
        if payload["solver"]["status"] != "INFEASIBLE":
            raise AssertionError(f"case {degree} is not certified infeasible")
        seen.add(degree)
        cases.append(
            {
                "root_negative_degree": degree,
                "source": str(path),
                "solver_status": payload["solver"]["status"],
                "solver_wall_time_seconds": payload["solver"][
                    "wall_time_seconds"
                ],
                "conflicts": payload["solver"]["conflicts"],
                "branches": payload["solver"]["branches"],
            }
        )
    if seen != set(range(6)):
        raise AssertionError((seen, set(range(6))))
    cases.sort(key=lambda case: case["root_negative_degree"])

    witness13_payload = json.loads(args.m13_witness.read_text())
    witness13 = matrix_from_payload(witness13_payload)
    profile13 = exact_profile(witness13)
    if len(witness13) != 13 or profile13["M"] != 20:
        raise AssertionError((len(witness13), profile13["M"]))
    witness14_payload = json.loads(args.m14_witness.read_text())
    witness14 = matrix_from_payload(witness14_payload)
    profile14 = exact_profile(witness14)
    if len(witness14) != 14 or profile14["M"] != 21:
        raise AssertionError((len(witness14), profile14["M"]))

    output = {
        "schema": "quadratic-signing-certified-m13-m14-v1",
        "classification": (
            "solver-certified exact finite values with explicit exhaustively "
            "verified upper witnesses; CP-SAT has no standalone proof object"
        ),
        "normalization": "M_n=max_x |sum_{i<j} a_ij x_i x_j|",
        "symmetry_argument": [
            "switching gauges all root edges to +1",
            "global signing negation followed by re-gauging complements all internal edges",
            "the same choice permits at most 33 negative internal edges",
            "a graph or its complement on 12 vertices has minimum degree at most 5",
            "permutation symmetry chooses a minimum-degree vertex as vertex 1",
            "the six root-degree cases 0 through 5 are disjoint and exhaustive",
        ],
        "cap18_cases": cases,
        "m13": {
            "certified_value": 20,
            "lower_chain": "six exhaustive cap-18 cases are infeasible; energies are even",
            "upper_source": str(args.m13_witness),
            "upper_matrix_sha256": stable_matrix_hash(witness13),
            "upper_profile": profile13,
        },
        "m14": {
            "certified_value": 21,
            "lower_chain": "M_14>=M_13=20 and all order-14 energies are odd",
            "upper_source": str(args.m14_witness),
            "upper_matrix_sha256": stable_matrix_hash(witness14),
            "upper_profile": profile14,
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print("certified M_13=20 and M_14=21")
    for case in cases:
        print(
            f"degree={case['root_negative_degree']} "
            f"status={case['solver_status']} "
            f"wall={case['solver_wall_time_seconds']:.6f}s"
        )
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
