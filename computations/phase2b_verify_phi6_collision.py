#!/usr/bin/env python3
"""Independently verify the order-10 equal-phi6/different-cap collision."""

from __future__ import annotations

import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

import numpy as np

from phase2_restriction_state_audit import (
    DisjointSet,
    class_map,
    root_gauge_code,
    signing_from_code,
)


ORDER = 10
CODES = (5850642905, 28771662001)
OUTPUT = Path("computations/results/phase2b_phi6_collision_n10.json")


def matrix_from_code(code: int) -> np.ndarray:
    matrix = np.ones((ORDER, ORDER), dtype=np.int64)
    np.fill_diagonal(matrix, 0)
    bit = 0
    for i in range(1, ORDER):
        for j in range(i + 1, ORDER):
            if code & (1 << bit):
                matrix[i, j] = matrix[j, i] = -1
            bit += 1
    return matrix


def exact_energy_profile(matrix: np.ndarray) -> dict[str, object]:
    codes = np.arange(1 << (len(matrix) - 1), dtype=np.uint64)
    bits = ((codes[:, None] >> np.arange(len(matrix) - 1, dtype=np.uint64)) & 1)
    spins = np.column_stack([np.ones(len(codes), dtype=np.int64), 1 - 2 * bits.astype(np.int64)])
    energies = np.einsum("bi,ij,bj->b", spins, matrix, spins) // 2
    histogram = Counter(map(int, energies))
    return {
        "minimum": int(energies.min()),
        "maximum": int(energies.max()),
        "cap": int(np.max(np.abs(energies))),
        "energy_histogram": {str(k): v for k, v in sorted(histogram.items())},
    }


def restriction_profile(matrix: np.ndarray, size: int,
                        labels: list[int]) -> list[int]:
    counts = [0] * (max(labels) + 1)
    for vertices in itertools.combinations(range(len(matrix)), size):
        child = matrix[np.ix_(vertices, vertices)]
        counts[labels[root_gauge_code(child)]] += 1
    return counts


def oriented_class_map(size: int) -> tuple[list[int], int]:
    """Switching/permutation classes without quotienting global negation."""
    bits = (size - 1) * (size - 2) // 2
    count = 1 << bits
    dsu = DisjointSet(count)
    for code in range(count):
        matrix = signing_from_code(code, size)
        for position in range(size - 1):
            permutation = list(range(size))
            permutation[position], permutation[position + 1] = (
                permutation[position + 1], permutation[position])
            moved = matrix[np.ix_(permutation, permutation)]
            dsu.union(code, root_gauge_code(moved))
    roots = [dsu.find(code) for code in range(count)]
    label_of = {root: label for label, root in enumerate(sorted(set(roots)))}
    return [label_of[root] for root in roots], len(label_of)


def positive_within_blowup(matrix: np.ndarray, multiplicity: int) -> np.ndarray:
    result = np.empty((len(matrix) * multiplicity,) * 2, dtype=np.int64)
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            block = result[
                i * multiplicity : (i + 1) * multiplicity,
                j * multiplicity : (j + 1) * multiplicity,
            ]
            if i == j:
                block.fill(1)
                np.fill_diagonal(block, 0)
            else:
                block.fill(matrix[i, j])
    return result


def bareiss_determinant(matrix: np.ndarray) -> int:
    a = [[int(value) for value in row] for row in matrix]
    n = len(a)
    sign = 1
    previous = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            swap = next((i for i in range(k + 1, n) if a[i][k]), None)
            if swap is None:
                return 0
            a[k], a[swap] = a[swap], a[k]
            sign *= -1
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                a[i][j] = (a[i][j] * pivot - a[i][k] * a[k][j]) // previous
        previous = pivot
    return sign * a[-1][-1]


def main() -> int:
    labels = {}
    class_counts = {}
    oriented_labels = {}
    oriented_class_counts = {}
    for size in (4, 5, 6):
        labels[size], class_counts[size] = class_map(size)
    for size in range(2, 7):
        oriented_labels[size], oriented_class_counts[size] = oriented_class_map(size)
    records = []
    for code in CODES:
        matrix = matrix_from_code(code)
        q = matrix + np.eye(ORDER, dtype=np.int64)
        shifted = q + 4 * np.eye(ORDER, dtype=np.int64)
        leading_minors = [
            bareiss_determinant(shifted[:k, :k]) for k in range(1, ORDER + 1)
        ]
        if not all(value > 0 for value in leading_minors):
            raise AssertionError("Q+4I failed the exact Sylvester certificate")
        records.append({
            "code": code,
            "matrix": matrix.tolist(),
            "matrix_sha256": hashlib.sha256(matrix.astype(np.int8).tobytes()).hexdigest(),
            "exact_energy_profile": exact_energy_profile(matrix),
            "restriction_class_histograms": {
                str(size): restriction_profile(matrix, size, labels[size])
                for size in (4, 5, 6)
            },
            "oriented_restriction_class_histograms": {
                str(size): restriction_profile(matrix, size, oriented_labels[size])
                for size in range(2, 7)
            },
            "Q_plus_4I_leading_principal_minors": leading_minors,
        })
    if records[0]["restriction_class_histograms"] != records[1]["restriction_class_histograms"]:
        raise AssertionError("phi6 profiles differ")
    if (records[0]["oriented_restriction_class_histograms"] !=
            records[1]["oriented_restriction_class_histograms"]):
        raise AssertionError("oriented profiles differ")
    if [r["exact_energy_profile"]["cap"] for r in records] != [19, 21]:
        raise AssertionError("cap collision mismatch")

    # Rational SDP-dual certificate for the first base Q=A+I.  Multiplying
    # diag(y) +/- Q by 100 gives the two integer matrices checked below.
    dual_numerators = np.asarray(
        [445, 490, 661, 668, 436, 645, 405, 427, 485, 513], dtype=np.int64
    )
    first_q = matrix_from_code(CODES[0]) + np.eye(ORDER, dtype=np.int64)
    dual_minors = {}
    for sign, label in ((-1, "diag_y_minus_Q"), (1, "diag_y_plus_Q")):
        integer_matrix = np.diag(dual_numerators) + sign * 100 * first_q
        minors = [
            bareiss_determinant(integer_matrix[:k, :k])
            for k in range(1, ORDER + 1)
        ]
        if not all(value > 0 for value in minors):
            raise AssertionError(f"Hadamard lift dual certificate failed: {label}")
        dual_minors[label] = minors
    # Independent finite audit of the blowup-preservation argument.
    blowups = [positive_within_blowup(matrix_from_code(code), 2) for code in CODES]
    blowup_profiles = []
    for blowup in blowups:
        blowup_profiles.append({
            str(size): restriction_profile(blowup, size, labels[size])
            for size in (4, 5, 6)
        })
    if blowup_profiles[0] != blowup_profiles[1]:
        raise AssertionError("twofold blowup phi6 profiles differ")

    output = {
        "schema": "quadratic-signing-phase2b-phi6-collision-v1",
        "classification": (
            "exact exhaustive spin caps, exact induced switching-class profiles, "
            "and exact positive-definiteness certificates"
        ),
        "order": ORDER,
        "search": {
            "classification": "deterministic nonexhaustive search; discovered pair independently verified exactly",
            "root_gauged_global_negation_representatives": 1 << 35,
            "distinct_codes_checked_before_collision": 555,
            "enumeration": "odd affine permutation modulo 2^35",
            "search_source": "computations/phase2_profile_collision_n8.cpp --sample 10 1000000",
        },
        "restriction_class_counts": {
            str(size): class_counts[size] for size in (4, 5, 6)
        },
        "oriented_restriction_class_counts": {
            str(size): oriented_class_counts[size] for size in range(2, 7)
        },
        "common_phi6_vector": [
            *records[0]["restriction_class_histograms"]["4"],
            *records[0]["restriction_class_histograms"]["5"],
            *records[0]["restriction_class_histograms"]["6"],
        ],
        "records": records,
        "balanced_all_positive_within_block_blowup": {
            "definition": (
                "replace every vertex by L twins, put +1 within each twin class, "
                "and retain the base sign between different classes"
            ),
            "valid_for_integer_L_at_least": 3,
            "first_cap_formula": "24*L^2-5*L",
            "second_cap_formula": "26*L^2-5*L",
            "cap_gap_formula": "2*L^2",
            "profile_statement": (
                "the two blowups have identical induced switching-class counts "
                "for every restriction order at most 6"
            ),
            "twofold_blowup_profile_independently_checked": blowup_profiles[0],
        },
        "sylvester_hadamard_low_scale_separation": {
            "orders": "k=4^r for every integer r>=1",
            "lift_definition": (
                "S_A(k)=A tensor H_k + I_10 tensor (H_k-diag(H_k)); "
                "define S_B(k) analogously"
            ),
            "profile_statement": (
                "S_A(k) and S_B(k) have identical oriented restriction profiles "
                "through order 6"
            ),
            "first_cap_upper": "cap(S_A(k)) <= (207/8)*k^(3/2)",
            "second_cap_lower": "cap(S_B(k)) >= 26*k^(3/2)",
            "cap_gap_lower": "cap(S_B(k))-cap(S_A(k)) >= k^(3/2)/8",
            "order_N": "N=10*k",
            "normalized_gap_lower": "1/(8*10^(3/2))",
            "dual_y_numerators_over_100": dual_numerators.tolist(),
            "dual_y_sum": "207/4",
            "dual_positive_definite_leading_minors_after_scaling_by_100":
                dual_minors,
            "positive_boolean_eigenvector_seed_for_H4": [-1, -1, -1, 1],
        },
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(f"verified phi6 collision caps 19,21; wrote {OUTPUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
