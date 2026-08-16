#!/usr/bin/env python3
"""Exact checks for the parity-half-cube replica hierarchy.

The analytic proof is Theorem 3.3 in ../theorems.md.  This script checks its
finite identities without relying on NetworkX:

* every proper selected codeword configuration has the same normalized
  column-pattern multiset in the two parity codes;
* the exact covering radii at r=3 and r=5;
* the alternating-binomial radius identity; and
* equality of the complete ambient t-point census through t=4, followed by
  separation at t=5, for the r=3 base pair.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from itertools import combinations, product
from math import comb
from pathlib import Path


def parity_coordinates(r: int, epsilon: int) -> list[int]:
    return [v for v in range(1 << r) if bin(v).count("1") % 2 == epsilon]


def parity_code(r: int, epsilon: int) -> tuple[list[int], int]:
    coordinates = parity_coordinates(r, epsilon)
    words = [0]
    for i in range(r):
        word = 0
        for column, v in enumerate(coordinates):
            if (v >> i) & 1:
                word |= 1 << column
        words.append(word)
    return words, len(coordinates)


def normalized_column_profile(
    words: list[int], length: int, labels: tuple[int, ...]
) -> Counter[tuple[int, ...]]:
    base = words[labels[0]]
    return Counter(
        tuple(((words[label] ^ base) >> column) & 1 for label in labels)
        for column in range(length)
    )


def verify_selected_profiles(r: int) -> int:
    code0, length = parity_code(r, 0)
    code1, length1 = parity_code(r, 1)
    assert length == length1
    checked = 0
    for size in range(1, r + 1):
        for labels in combinations(range(r + 1), size):
            assert normalized_column_profile(code0, length, labels) == (
                normalized_column_profile(code1, length, labels)
            )
            checked += 1
    return checked


def covering_radius(code: list[int], length: int) -> int:
    return max(
        min(bin(root ^ word).count("1") for word in code)
        for root in range(1 << length)
    )


def s_value(r: int, epsilon: int) -> int:
    order = r + 1
    return sum(
        comb(r, weight) * max(weight, order - weight)
        for weight in range(r + 1)
        if weight % 2 == epsilon
    )


def t_census(code: set[int], length: int, replicas: int) -> Counter[tuple]:
    cube = range(1 << length)
    answer: Counter[tuple] = Counter()
    for states in product(cube, repeat=replicas):
        membership = tuple(int(state in code) for state in states)
        distances = tuple(
            bin(states[i] ^ states[j]).count("1")
            for i in range(replicas)
            for j in range(i + 1, replicas)
        )
        answer[membership + distances] += 1
    return answer


def verify_r3_censuses() -> list[dict]:
    code0, length = parity_code(3, 0)
    code1, _ = parity_code(3, 1)
    results = []
    for replicas in range(1, 6):
        census0 = t_census(set(code0), length, replicas)
        census1 = t_census(set(code1), length, replicas)
        equal = census0 == census1
        assert equal == (replicas <= 4)
        results.append(
            {
                "replicas": replicas,
                "equal": equal,
                "cells_parity_0": len(census0),
                "cells_parity_1": len(census1),
                "total_tuples": (1 << length) ** replicas,
            }
        )
    return results


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "extremal_information/experiments/"
            "code_replica_hierarchy_results.json"
        ),
    )
    args = parser.parse_args()

    profile_checks = []
    for r in (3, 5, 7, 9):
        profile_checks.append(
            {
                "r": r,
                "order": r + 1,
                "block_length": 1 << (r - 1),
                "proper_label_subsets_checked": verify_selected_profiles(r),
            }
        )

    exact_radii = []
    for r in (3, 5):
        row = {"r": r, "block_length": 1 << (r - 1)}
        radii = []
        for epsilon in (0, 1):
            code, length = parity_code(r, epsilon)
            radii.append(covering_radius(code, length))
            row[f"S_{epsilon}"] = s_value(r, epsilon)
        p = ((r - 1) // 2) % 2
        assert s_value(r, p) - s_value(r, 1 - p) == comb(r - 1, (r - 1) // 2)
        assert radii[p] == s_value(r, p) // (r + 1)
        assert radii[p] > radii[1 - p]
        row["covering_radii"] = radii
        row["larger_radius_parity"] = p
        row["normalized_cartesian_power_gap"] = abs(radii[0] - radii[1]) / length
        exact_radii.append(row)

    result = {
        "schema": "extremal-information-code-replica-hierarchy-v1",
        "selected_profile_checks": profile_checks,
        "exact_radius_checks": exact_radii,
        "r3_complete_ambient_censuses": verify_r3_censuses(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
