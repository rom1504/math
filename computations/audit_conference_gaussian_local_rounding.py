#!/usr/bin/env python3
"""Exact audit of Gaussian top-space rounding and its local-search barrier."""

from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from two_fiber_cyclic_conference import CERTIFICATES, circulant, verify


def conference_matrix(k: int) -> np.ndarray:
    data = CERTIFICATES[k]
    verify(k, data["a"], data["c"])
    matrix_a = circulant(data["a"])
    matrix_c = circulant(data["c"])
    return np.block([[matrix_a, matrix_c], [matrix_c.T, -matrix_a]])


def root_trap(matrix: np.ndarray, root: int) -> dict[str, object]:
    order = len(matrix)
    q = order - 1
    spin = matrix[:, root].copy()
    spin[root] = 1
    image = matrix @ spin
    fields = spin * image
    if int(spin @ image // 2) != q:
        raise AssertionError("root trap energy failed")
    if Counter(map(int, fields)) != Counter({1: q, q: 1}):
        raise AssertionError("root trap local fields failed")
    # P_+ spin=(spin+S spin/sqrt(q))/2.  Its coordinate signs equal spin
    # exactly iff q+sqrt(q)>0 at the root and 1+sqrt(q)>0 elsewhere.
    if np.any(fields <= 0):
        raise AssertionError("root trap is not a strict local maximum")
    if not np.array_equal(np.sign(image).astype(np.int64), spin):
        raise AssertionError("synchronous local update does not fix trap")
    return {
        "root": root,
        "positive_bit_hex_little_endian": format(
            sum(1 << index for index, value in enumerate(spin) if value > 0), "x"
        ),
        "exact_energy": q,
        "normalized_energy": q / order**1.5,
        "local_field_histogram": {"1": q, str(q): 1},
        "strict_one_flip_local_maximum": True,
        "fixed_by_synchronous_sign_Sx": True,
        "is_sign_of_top_projection_P_plus_x": True,
    }


def flip_set_audit(matrix: np.ndarray, root: int = 0) -> dict[str, object]:
    order = len(matrix)
    q = order - 1
    trap = matrix[:, root].copy()
    trap[root] = 1
    gauged = trap[:, None] * matrix * trap[None, :]
    other = np.asarray([index for index in range(order) if index != root])
    core = gauged[np.ix_(other, other)]
    if not np.array_equal(core @ np.ones(q, dtype=np.int64), np.zeros(q, dtype=np.int64)):
        raise AssertionError("normalized core is not row-balanced")
    if not np.array_equal(core @ core, q * np.eye(q, dtype=np.int64) - np.ones((q, q), dtype=np.int64)):
        raise AssertionError("normalized core square identity failed")
    positive = (core > 0).astype(np.int64)
    np.fill_diagonal(positive, 0)
    degree = (q - 1) // 2
    adjacent_common = (q - 5) // 4
    nonadjacent_common = (q - 1) // 4
    if not np.all(positive.sum(axis=1) == degree):
        raise AssertionError("positive core graph degree failed")
    square = positive @ positive
    for left in range(q):
        for right in range(left + 1, q):
            expected = adjacent_common if positive[left, right] else nonadjacent_common
            if square[left, right] != expected:
                raise AssertionError("conference graph common-neighbor count failed")
    triangle = next(
        (
            vertices
            for vertices in itertools.combinations(range(q), 3)
            if all(positive[i, j] for i, j in itertools.combinations(vertices, 2))
        ),
        None,
    )
    baseline = q
    result: dict[str, object] = {
        "root": root,
        "core_order": q,
        "positive_core_graph_parameters": {
            "v": q,
            "degree": degree,
            "lambda": adjacent_common,
            "mu": nonadjacent_common,
        },
        "exact_improvement_criterion": (
            "a core flip set F of size m improves iff its induced positive-edge "
            "count p(F) satisfies p(F)>m^2/4"
        ),
        "exact_local_fields_for_general_core_flip_set": (
            "t_root=q-2m; t_v=2m+1-4d_F(v) outside F; "
            "t_v=4d_F(v)-2m+1 inside F"
        ),
        "proved_terminal_basin_radius": (
            "for q>=9, every nonempty coordinate local maximum has "
            "m>=(3+sqrt(q))/2"
        ),
        "one_flip_strict_improvement": False,
        "two_flip_strict_improvement": False,
    }
    if triangle is None:
        if q != 5:
            raise AssertionError("triangle unexpectedly absent")
        result.update({
            "smallest_strictly_improving_flip_set_size": None,
            "classification": (
                "positive core is the triangle-free conference graph on five "
                "vertices; Mantel excludes every improving flip set"
            ),
        })
    else:
        full_vertices = [int(other[index]) for index in triangle]
        candidate = np.ones(order, dtype=np.int64)
        candidate[full_vertices] = -1
        energy = int(candidate @ gauged @ candidate // 2)
        if energy != baseline + 6:
            raise AssertionError("triangle flip energy failed")
        fields = candidate * (gauged @ candidate)
        outside = [index for index in range(q) if index not in triangle]
        neighbor_counts = Counter(
            int(sum(positive[index, vertex] for vertex in triangle))
            for index in outside
        )
        rho = neighbor_counts.get(3, 0)
        expected_neighbor_counts = {
            0: (q - 9) // 4 - rho,
            1: 6 + 3 * rho,
            2: 3 * (q - 9) // 4 - 3 * rho,
            3: rho,
        }
        expected_counter = Counter(
            {key: value for key, value in expected_neighbor_counts.items() if value}
        )
        if neighbor_counts != expected_counter:
            raise AssertionError((neighbor_counts, expected_neighbor_counts))
        negative_count = int(np.sum(fields < 0))
        if negative_count != 3 * (q - 9) // 4 - 2 * rho:
            raise AssertionError("post-triangle negative count failed")
        result.update({
            "smallest_strictly_improving_flip_set_size": 3,
            "positive_triangle_core_indices": list(triangle),
            "positive_triangle_full_indices": full_vertices,
            "energy_before": baseline,
            "energy_after_triangle_flip": energy,
            "exact_energy_gain": 6,
            "classification": "every normalized conference core of order at least nine has a positive triangle",
            "post_triangle_local_field_histogram": {
                str(key): value
                for key, value in sorted(Counter(map(int, fields)).items())
            },
            "triangle_K4_extension_count_rho": rho,
            "outside_triangle_positive_neighbor_histogram": {
                str(key): value for key, value in sorted(neighbor_counts.items())
            },
            "post_triangle_strictly_improving_one_flip_count": negative_count,
            "universal_negative_count_lower_bound": (q - 9) // 4,
            "best_immediate_one_flip_gain": 10 if rho else 2,
        })
    return result


def exhaustive_local_maximum_audit(matrix: np.ndarray) -> dict[str, object]:
    order = len(matrix)
    codes = np.arange(1 << (order - 1), dtype=np.uint64)
    positions = np.arange(order - 1, dtype=np.uint64)
    tails = 1 - 2 * ((codes[:, None] >> positions) & 1).astype(np.int8)
    spins = np.column_stack((np.ones(len(codes), dtype=np.int8), tails))
    images = spins @ matrix
    fields = spins * images
    energies = np.sum(fields, axis=1) // 2
    positive_local = np.all(fields > 0, axis=1)
    local_energies = energies[positive_local]
    if len(local_energies) == 0 or int(np.min(local_energies)) != order - 1:
        raise AssertionError("sharp local-minimum audit failed")
    return {
        "order": order,
        "projective_spin_count": len(spins),
        "strict_positive_local_maximum_count": int(np.sum(positive_local)),
        "minimum_positive_local_maximum_energy": int(np.min(local_energies)),
        "energy_histogram": {
            str(key): value
            for key, value in sorted(Counter(map(int, local_energies)).items())
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    records = []
    for k in (1, 2, 3):
        matrix = conference_matrix(k)
        order = len(matrix)
        q = order - 1
        expected_energy = order * q / math.pi * math.asin(1 / math.sqrt(q))
        record = {
            "k": k,
            "order": order,
            "conference_multiplier": q,
            "exact_gaussian_rounding_expectation_formula": (
                "N*(N-1)/pi*asin(1/sqrt(N-1))"
            ),
            "gaussian_rounding_expected_energy": expected_energy,
            "gaussian_rounding_expected_normalized_energy": expected_energy / order**1.5,
            "root_trap": root_trap(matrix, 0),
            "flip_set_audit": flip_set_audit(matrix, 0),
        }
        if order <= 18:
            record["exhaustive_local_maximum_audit"] = exhaustive_local_maximum_audit(
                matrix
            )
        records.append(record)

    payload = {
        "schema": "quadratic-signing-conference-gaussian-local-rounding-v1",
        "classification": (
            "proved exact Gaussian expectation and sharp universal pointwise "
            "barrier for one-flip local improvement"
        ),
        "universal_results": {
            "gaussian_expectation": (
                "E H(sign(g))=N*(N-1)/pi*asin(1/sqrt(N-1)) "
                "for g~N(0,P_plus)"
            ),
            "asymptotic_expected_constant": "1/pi",
            "sharp_local_maximum_lower_bound": (
                "every positive one-flip local maximum has H>=N-1"
            ),
            "sharp_root_trap": (
                "x=e_r+S e_r has H=N-1, local fields {N-1,1,...,1}, "
                "x=sign(P_plus x), and sign(Sx)=x"
            ),
            "post_triangle_escape": (
                "a positive-triangle flip gains 6 and creates at least "
                "(N-10)/4 negative local fields"
            ),
            "bounded_radius_improving_supply": (
                "for every fixed core flip-set size m>=3 there are Omega(N) "
                "improving additions; the proved uniform lower bound is "
                "Omega(N/m) while m=o(sqrt(N))"
            ),
            "terminal_basin_radius": (
                "for N>=10, every nonempty normalized terminal flip set has "
                "m>=(3+sqrt(N-1))/2"
            ),
        },
        "finite_exact_audits": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps(payload, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
