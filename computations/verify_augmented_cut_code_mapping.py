#!/usr/bin/env python3
"""Verify the exact augmented-cut-code mapping on certified minimizers.

Let Cut(K_n) be the binary cut code and let C_n^+ add the all-one edge word.
For a signing word a and its quadratic cap Q(a),

    Q(a) = binom(n,2) - 2 d(a,C_n^+).

Consequently M_n=binom(n,2)-2 rho(C_n^+).  The script verifies the fixed-word
identity by explicit codeword enumeration and records the covering radii that
follow from the independently certified exact values through order 14.
"""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np

from exact_mn_milp import exact_profile, stable_matrix_hash


DEFAULT_WITNESSES = {
    3: ("computations/results/exact_m3.json", "matrix", 3),
    4: ("computations/results/exact_m4.json", "matrix", 4),
    5: ("computations/results/exact_m5.json", "matrix", 4),
    6: ("computations/results/exact_m6.json", "matrix", 5),
    7: ("computations/results/exact_m7.json", "matrix", 9),
    8: ("computations/results/exact_m8.json", "matrix", 10),
    9: ("computations/results/exact_m9.json", "matrix", 12),
    10: ("computations/results/exact_m10.json", "matrix", 13),
    11: ("computations/results/heuristic_m11.json", "matrix", 17),
    12: (
        "computations/results/extension_nested_m11_to_12.json",
        "parent_matrix",
        18,
    ),
    13: (
        "computations/results/bridge_6_7_sign1_cap20.json",
        "parent_matrix",
        20,
    ),
    14: (
        "computations/results/conference_completion_m13.json",
        "conference_matrix",
        21,
    ),
}


def distance_to_augmented_cut_code(matrix: np.ndarray) -> tuple[int, int]:
    n = len(matrix)
    edges = tuple(itertools.combinations(range(n), 2))
    word = np.asarray([matrix[i, j] == -1 for i, j in edges], dtype=np.int8)
    minimum = len(edges) + 1
    codewords = set()
    for mask in range(1 << (n - 1)):
        labels = np.zeros(n, dtype=np.int8)
        labels[1:] = [(mask >> bit) & 1 for bit in range(n - 1)]
        cut = np.asarray([labels[i] ^ labels[j] for i, j in edges], dtype=np.int8)
        for constant in (0, 1):
            codeword = cut ^ constant
            codewords.add(codeword.tobytes())
            minimum = min(minimum, int(np.count_nonzero(word ^ codeword)))
    if len(codewords) != 1 << n:
        raise AssertionError((n, len(codewords), 1 << n))
    return minimum, len(codewords)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    records = []
    for n, (source_text, key, exact_value) in DEFAULT_WITNESSES.items():
        source = Path(source_text)
        payload = json.loads(source.read_text())
        matrix = np.asarray(payload[key], dtype=np.int8)
        profile = exact_profile(matrix)
        if profile["M"] != exact_value:
            raise AssertionError((n, profile["M"], exact_value))
        distance, code_size = distance_to_augmented_cut_code(matrix)
        edge_count = n * (n - 1) // 2
        if edge_count - 2 * distance != exact_value:
            raise AssertionError((n, edge_count, distance, exact_value))
        records.append(
            {
                "n": n,
                "edge_count": edge_count,
                "augmented_cut_code_dimension": n,
                "augmented_cut_code_size": code_size,
                "certified_M_n": exact_value,
                "covering_radius": distance,
                "source": source_text,
                "matrix_key": key,
                "matrix_sha256": stable_matrix_hash(matrix),
            }
        )
        print(
            f"n={n} M={exact_value} rho(C+)=({edge_count}-{exact_value})/2="
            f"{distance} code_size={code_size}"
        )
    output = {
        "schema": "quadratic-signing-augmented-cut-covering-radius-v1",
        "classification": (
            "proved algebraic mapping, explicit fixed-word verification, and "
            "arithmetic consequences of solver-certified exact M_n values"
        ),
        "mapping": "M_n=binom(n,2)-2*rho(Cut(K_n)+span{all-one edge word})",
        "dual_description": (
            "the dual is the even-weight subcode of the binary cycle space of K_n"
        ),
        "records": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
