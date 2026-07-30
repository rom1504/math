#!/usr/bin/env python3
"""Construct the GF(9) Paley conference matrix and audit its restrictions."""

from __future__ import annotations

import argparse
import itertools
import json
from collections import Counter
from pathlib import Path

import numpy as np

from analyze_equal_split_partitions import add_orbit_class, classify
from exact_mn_milp import exact_profile, stable_matrix_hash


# Represent a+b*t by (a,b), with t^2=-1=2 in F_3.
ELEMENTS = tuple((a, b) for a in range(3) for b in range(3))


def add(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    return ((x[0] + y[0]) % 3, (x[1] + y[1]) % 3)


def negate(x: tuple[int, int]) -> tuple[int, int]:
    return ((-x[0]) % 3, (-x[1]) % 3)


def multiply(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
    return (
        (x[0] * y[0] + 2 * x[1] * y[1]) % 3,
        (x[0] * y[1] + x[1] * y[0]) % 3,
    )


def power(x: tuple[int, int], exponent: int) -> tuple[int, int]:
    result = (1, 0)
    base = x
    while exponent:
        if exponent & 1:
            result = multiply(result, base)
        base = multiply(base, base)
        exponent //= 2
    return result


def character(x: tuple[int, int]) -> int:
    if x == (0, 0):
        return 0
    fourth = power(x, 4)
    if fourth == (1, 0):
        return 1
    if fourth == (2, 0):
        return -1
    raise AssertionError((x, fourth))


def conference() -> np.ndarray:
    matrix = np.zeros((10, 10), dtype=np.int8)
    matrix[0, 1:] = matrix[1:, 0] = 1
    for i, x in enumerate(ELEMENTS, start=1):
        for j, y in enumerate(ELEMENTS, start=1):
            if i != j:
                matrix[i, j] = character(add(x, negate(y)))
    if not np.array_equal(
        matrix.astype(np.int64) @ matrix.astype(np.int64),
        9 * np.eye(10, dtype=np.int64),
    ):
        raise AssertionError("conference identity failed")
    return matrix


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("classes", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    matrix = conference()
    profile = exact_profile(matrix)
    class_payload = json.loads(args.classes.read_text())
    classes: list[dict[str, set[bytes]]] = []
    for row in class_payload["classes"]:
        add_orbit_class(
            classes, np.asarray(row["representative_matrix"], dtype=np.int8)
        )

    deletion9_caps: Counter[int] = Counter()
    for deleted in range(10):
        keep = [vertex for vertex in range(10) if vertex != deleted]
        deletion9_caps[exact_profile(matrix[np.ix_(keep, keep)])["M"]] += 1
    deletion8_caps: Counter[int] = Counter()
    deletion8_classes: Counter[str] = Counter()
    examples = {}
    for deleted in itertools.combinations(range(10), 2):
        keep = [vertex for vertex in range(10) if vertex not in deleted]
        principal = matrix[np.ix_(keep, keep)]
        cap = exact_profile(principal)["M"]
        deletion8_caps[cap] += 1
        if cap == 10:
            class_index, sign = classify(classes, principal)
            key = f"class{class_index}{sign}"
            deletion8_classes[key] += 1
            examples.setdefault(
                key,
                {
                    "deleted_vertices": list(deleted),
                    "matrix": [[int(value) for value in row] for row in principal],
                    "matrix_sha256": stable_matrix_hash(principal),
                },
            )
    output = {
        "schema": "quadratic-signing-conference-order10-gf9-v1",
        "classification": (
            "proved GF(9) construction, exact integer identities, and exhaustive finite profiles"
        ),
        "field_model": "GF(9)=GF(3)[t]/(t^2+1)",
        "conference_matrix": [[int(value) for value in row] for row in matrix],
        "conference_matrix_sha256": stable_matrix_hash(matrix),
        "conference_identity": "C^2=9I",
        "conference_profile": profile,
        "order9_deletion_cap_counts": {
            str(cap): count for cap, count in sorted(deletion9_caps.items())
        },
        "order8_double_deletion_cap_counts": {
            str(cap): count for cap, count in sorted(deletion8_caps.items())
        },
        "order8_minimizer_class_counts": dict(sorted(deletion8_classes.items())),
        "order8_class_examples": examples,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"C10 cap={profile['M']} deletions9={dict(deletion9_caps)} "
        f"deletions8={dict(deletion8_caps)} classes={dict(deletion8_classes)}"
    )
    print(f"hash={output['conference_matrix_sha256']}")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
