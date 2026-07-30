#!/usr/bin/env python3
"""Structural and cross-order statistics for saved signing witnesses.

All norms and restriction statistics are evaluated exhaustively.  Conclusions
about a saved witness are exact; conclusions about all minimizers require a
separate enumeration and are not inferred by this script.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from itertools import combinations
from pathlib import Path

import numpy as np

from exact_mn_milp import exact_profile, projective_spins, stable_matrix_hash


KNOWN_EXACT_M = {
    3: 3,
    4: 4,
    5: 4,
    6: 5,
    7: 9,
    8: 10,
    9: 12,
    10: 13,
    11: 17,
    12: 18,
}


def boolean_cap(matrix: np.ndarray) -> int:
    spins = projective_spins(len(matrix)).astype(np.int64)
    energies = np.einsum("bi,ij,bj->b", spins, matrix.astype(np.int64), spins) // 2
    return int(np.abs(energies).max())


def restriction_statistics(matrix: np.ndarray) -> dict[str, object]:
    n = len(matrix)
    by_size: dict[str, object] = {}
    for k in range(3, n):
        values = []
        for subset in combinations(range(n), k):
            child = matrix[np.ix_(subset, subset)]
            values.append(boolean_cap(child))
        counts = Counter(values)
        entry: dict[str, object] = {
            "minimum": min(values),
            "maximum": max(values),
            "histogram": {str(v): counts[v] for v in sorted(counts)},
        }
        if k in KNOWN_EXACT_M:
            entry["known_M"] = KNOWN_EXACT_M[k]
            entry["optimal_restriction_count"] = counts[KNOWN_EXACT_M[k]]
        by_size[str(k)] = entry
    return by_size


def partition_statistics(matrix: np.ndarray) -> dict[str, object]:
    n = len(matrix)
    output: dict[str, object] = {}
    for k in range(3, n // 2 + 1):
        other = n - k
        if other < 3:
            continue
        records = []
        for subset in combinations(range(n), k):
            complement = tuple(i for i in range(n) if i not in subset)
            left = boolean_cap(matrix[np.ix_(subset, subset)])
            right = boolean_cap(matrix[np.ix_(complement, complement)])
            records.append((left, right))
        pair_counts = Counter(records)
        entry: dict[str, object] = {
            "partition_count": len(records),
            "child_cap_pairs": {
                f"{a},{b}": pair_counts[(a, b)] for a, b in sorted(pair_counts)
            },
            "minimum_cap_sum": min(a + b for a, b in records),
            "minimum_two_thirds_sum": min(
                a ** (2.0 / 3.0) + b ** (2.0 / 3.0) for a, b in records
            ),
        }
        if k in KNOWN_EXACT_M and other in KNOWN_EXACT_M:
            target = (KNOWN_EXACT_M[k], KNOWN_EXACT_M[other])
            entry["tight_optimal_child_partition_count"] = pair_counts[target]
        output[f"{k}+{other}"] = entry
    return output


def analyze(path: Path) -> dict[str, object]:
    payload = json.loads(path.read_text())
    key = (
        "matrix"
        if "matrix" in payload
        else ("parent_matrix" if "parent_matrix" in payload else "conference_matrix")
    )
    matrix = np.asarray(payload[key], dtype=np.int8)
    profile = exact_profile(matrix)
    expected_hash = payload.get(
        "matrix_sha256",
        payload.get("parent_matrix_sha256", payload.get("conference_matrix_sha256")),
    )
    if stable_matrix_hash(matrix) != expected_hash:
        raise AssertionError(f"matrix hash mismatch in {path}")
    n = len(matrix)
    square = matrix.astype(np.int64) @ matrix.astype(np.int64)
    spectral_defect = square - (n - 1) * np.eye(n, dtype=np.int64)
    eigenvalues = np.linalg.eigvalsh(matrix.astype(float))
    return {
        "source": str(path),
        "classification": payload.get("classification", "solver-certified exact optimization"),
        "n": n,
        "M": profile["M"],
        "normalized_M": profile["M"] / (n ** 1.5),
        "P": profile["P"],
        "Q": profile["Q"],
        "absolute_projective_ground_count": (
            profile["projective_top_count"]
            if profile["P"] == profile["M"]
            else 0
        )
        + (
            profile["projective_bottom_count"]
            if profile["Q"] == profile["M"]
            else 0
        ),
        "row_sums_sorted": sorted(profile["row_sums"]),
        "row_square_sum": int(np.sum(np.asarray(profile["row_sums"], dtype=np.int64) ** 2)),
        "operator_norm": float(np.max(np.abs(eigenvalues))),
        "eigenvalues": [float(v) for v in eigenvalues],
        "trace_A4": int(np.trace(square @ square)),
        "conference_defect_frobenius_squared": int(np.sum(spectral_defect**2)),
        "restrictions": restriction_statistics(matrix),
        "partitions": partition_statistics(matrix),
    }


def composition_defects(exact_values: dict[int, int]) -> list[dict[str, object]]:
    records = []
    for total in sorted(exact_values):
        for m in sorted(exact_values):
            n = total - m
            if m > n or n not in exact_values:
                continue
            defect = exact_values[total] ** (2.0 / 3.0) - exact_values[m] ** (
                2.0 / 3.0
            ) - exact_values[n] ** (2.0 / 3.0)
            records.append(
                {
                    "m": m,
                    "n": n,
                    "total": total,
                    "defect": defect,
                    "defect_over_total": defect / total,
                }
            )
    return records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("inputs", nargs="+", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    analyses = [analyze(path) for path in args.inputs]
    exact_values = dict(KNOWN_EXACT_M)
    for record in analyses:
        if record["n"] not in exact_values and "heuristic" not in record["classification"]:
            exact_values[int(record["n"])] = int(record["M"])
    payload = {
        "schema": "quadratic-signing-family-analysis-v1",
        "classification": (
            "exact exhaustive statistics for the saved witnesses; "
            "not an enumeration of all minimizers"
        ),
        "known_exact_M": {str(k): v for k, v in sorted(exact_values.items())},
        "witnesses": analyses,
        "two_thirds_composition_defects": composition_defects(exact_values),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    for record in analyses:
        deletion = record["restrictions"].get(str(record["n"] - 1), {})
        print(
            f"n={record['n']} M={record['M']} normalized={record['normalized_M']:.9f} "
            f"grounds={record['absolute_projective_ground_count']} "
            f"delete_range={deletion.get('minimum')}..{deletion.get('maximum')} "
            f"optimal_deletions={deletion.get('optimal_restriction_count')}"
        )
    positive = [r for r in payload["two_thirds_composition_defects"] if r["defect"] > 0]
    print(f"positive two-thirds defects: {len(positive)}")
    for record in positive:
        print(
            f"  {record['m']}+{record['n']}: defect={record['defect']:+.12f}"
        )
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
