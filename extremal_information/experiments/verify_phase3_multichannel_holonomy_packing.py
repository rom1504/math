#!/usr/bin/env python3
"""Finite falsification checks for the multichannel holonomy packing.

The asymptotic existence of the host code in MP.1 is analytic.  This script
checks three finite claims used by MP.2:

1. the two-fragment Cayley distance equals the subset formula (MP.2);
2. the response/Hausdorff comparison (MP.3a) on exhaustive small families;
3. all 2-subspaces of the Reed--Muller [16,5,8] code give distinct profiles
   with the predicted witness separation at least d-2k=4.

Only the standard library is used.  The JSON output is deterministic.
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
RESULT_PATH = HERE / "phase3_multichannel_holonomy_packing_results.json"


def popcount(value: int) -> int:
    # The project environment intentionally supports the older system Python.
    return bin(value).count("1")


def span(basis: tuple[int, ...]) -> frozenset[int]:
    values = {0}
    for vector in basis:
        values |= {x ^ vector for x in tuple(values)}
    return frozenset(values)


def profile_value(basis: tuple[int, ...], u: int) -> int:
    best = popcount(u)
    k = len(basis)
    for mask in range(1 << k):
        word = 0
        used = 0
        for j, vector in enumerate(basis):
            if (mask >> j) & 1:
                word ^= vector
                used += 1
        best = min(best, 2 * used + popcount(u ^ word))
    return best


def distance_to_code(u: int, code: frozenset[int]) -> int:
    return min(popcount(u ^ c) for c in code)


def hausdorff_distance(code: frozenset[int], other: frozenset[int]) -> int:
    return max(
        max(distance_to_code(c, other) for c in code),
        max(distance_to_code(c, code) for c in other),
    )


def subset_word_distances(columns: tuple[int, ...]) -> dict[int, int]:
    """Exact word distance for a small binary Cayley generating family."""

    distances: dict[int, int] = {0: 0}
    for column in columns:
        old = tuple(distances.items())
        for syndrome, weight in old:
            target = syndrome ^ column
            candidate = weight + 1
            if candidate < distances.get(target, candidate + 1):
                distances[target] = candidate
    return distances


def exhaustive_formula_check() -> dict[str, int]:
    # D=5, k=2 is large enough to have genuinely interacting channels but
    # small enough to check every ordered independent basis.
    D = 5
    k = 2
    checked_bases = 0
    checked_queries = 0
    max_metric_slack = 0
    bases: list[tuple[int, int]] = []

    for v1 in range(1, 1 << D):
        for v2 in range(1, 1 << D):
            if v2 == v1:
                continue
            basis = (v1, v2)
            bases.append(basis)

            # B, the zero lifts, and the V lifts in W direct-sum Q.
            kernel = tuple(1 << i for i in range(D))
            zero_lifts = tuple(1 << (D + j) for j in range(k))
            moved_lifts = tuple(v | (1 << (D + j)) for j, v in enumerate(basis))
            distances = subset_word_distances(kernel + zero_lifts + moved_lifts)

            code = span(basis)
            for u in range(1 << D):
                expected = profile_value(basis, u)
                actual = distances[u]
                if expected != actual:
                    raise AssertionError((basis, u, expected, actual))
                d_code = distance_to_code(u, code)
                if not d_code <= expected <= d_code + 2 * k:
                    raise AssertionError((basis, u, d_code, expected))
                checked_queries += 1
            checked_bases += 1

    # Exhaustively check MP.3a on one canonical basis for every 2-subspace of
    # F_2^5 (155 subspaces, 11,935 pairs).
    subspaces: dict[frozenset[int], tuple[int, int]] = {}
    for basis in bases:
        code = span(basis)
        subspaces.setdefault(code, basis)

    canonical = list(subspaces.items())
    metric_pairs = 0
    for index, (code, basis) in enumerate(canonical):
        values = [profile_value(basis, u) for u in range(1 << D)]
        for other, other_basis in canonical[index + 1 :]:
            other_values = [profile_value(other_basis, u) for u in range(1 << D)]
            response_distance = max(abs(a - b) for a, b in zip(values, other_values))
            hdist = hausdorff_distance(code, other)
            slack = abs(response_distance - hdist)
            if slack > 2 * k:
                raise AssertionError((basis, other_basis, response_distance, hdist))
            max_metric_slack = max(max_metric_slack, slack)
            metric_pairs += 1

    return {
        "ordered_independent_bases": checked_bases,
        "formula_queries": checked_queries,
        "grassmannian_subspaces": len(canonical),
        "grassmannian_pairs": metric_pairs,
        "maximum_observed_metric_slack": max_metric_slack,
    }


def rm_1_4_generators() -> tuple[int, ...]:
    """Truth tables of 1,x_0,x_1,x_2,x_3 on F_2^4."""

    generators = [(1 << 16) - 1]
    for coordinate in range(4):
        mask = 0
        for point in range(16):
            if (point >> coordinate) & 1:
                mask |= 1 << point
        generators.append(mask)
    return tuple(generators)


def linear_combination(generators: tuple[int, ...], coefficient: int) -> int:
    value = 0
    for i, generator in enumerate(generators):
        if (coefficient >> i) & 1:
            value ^= generator
    return value


def rm_host_packing_check() -> dict[str, int]:
    D = 16
    k = 2
    generators = rm_1_4_generators()
    host = span(generators)
    minimum_distance = min(popcount(x) for x in host if x)
    if len(host) != 32 or minimum_distance != 8:
        raise AssertionError((len(host), minimum_distance))

    # Enumerate coefficient-space 2-subspaces, retaining one ordered basis.
    coefficient_subspaces: dict[frozenset[int], tuple[int, int]] = {}
    for a, b in itertools.combinations(range(1, 1 << len(generators)), 2):
        coefficient_subspaces.setdefault(span((a, b)), (a, b))
    if len(coefficient_subspaces) != 155:  # Gaussian binomial [5 choose 2]_2.
        raise AssertionError(len(coefficient_subspaces))

    records: list[tuple[frozenset[int], tuple[int, int], frozenset[int]]] = []
    for coefficient_code, coefficient_basis in coefficient_subspaces.items():
        basis = tuple(linear_combination(generators, c) for c in coefficient_basis)
        records.append((coefficient_code, basis, span(basis)))

    pair_checks = 0
    minimum_witness_gap = D
    minimum_cross_distance = D
    for index, (coefficient_code, basis, code) in enumerate(records):
        for other_coefficients, other_basis, other_code in records[index + 1 :]:
            witness_coefficient = min(coefficient_code - other_coefficients)
            witness = linear_combination(generators, witness_coefficient)
            cross_distance = distance_to_code(witness, other_code)
            if cross_distance < minimum_distance:
                raise AssertionError((coefficient_code, other_coefficients, cross_distance))
            gap = profile_value(other_basis, witness) - profile_value(basis, witness)
            if gap < minimum_distance - 2 * k:
                raise AssertionError((basis, other_basis, cross_distance, gap))
            minimum_witness_gap = min(minimum_witness_gap, gap)
            minimum_cross_distance = min(minimum_cross_distance, cross_distance)
            pair_checks += 1

    return {
        "ambient_length": D,
        "host_dimension": len(generators),
        "host_minimum_distance": minimum_distance,
        "two_subspaces": len(records),
        "subspace_pairs": pair_checks,
        "minimum_cross_distance": minimum_cross_distance,
        "minimum_profile_witness_gap": minimum_witness_gap,
        "theoretical_gap_d_minus_2k": minimum_distance - 2 * k,
    }


def main() -> None:
    result = {
        "status": "passed",
        "exhaustive_small": exhaustive_formula_check(),
        "reed_muller_host": rm_host_packing_check(),
    }
    RESULT_PATH.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
