#!/usr/bin/env python3
"""Finite checks for the carrier-capacity law and rank-metric example."""

from __future__ import annotations

import itertools
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
RESULT_PATH = HERE / "phase3_carrier_capacity_results.json"


def span(basis: tuple[int, ...]) -> frozenset[int]:
    values = {0}
    for vector in basis:
        values |= {x ^ vector for x in tuple(values)}
    return frozenset(values)


def cycle_distance(x: int, y: int, size: int) -> int:
    delta = abs(x - y)
    return min(delta, size - delta)


def finite_metric_check() -> dict[str, int]:
    """Exhaust CC.3 over all {0,1}-cost carriers in the five-cycle."""

    size = 5
    weighted_carriers: list[tuple[frozenset[int], dict[int, int], tuple[int, ...]]] = []
    for carrier_mask in range(1, 1 << size):
        carrier = frozenset(i for i in range(size) if (carrier_mask >> i) & 1)
        carrier_list = sorted(carrier)
        for cost_mask in range(1 << len(carrier_list)):
            costs = {c: (cost_mask >> j) & 1 for j, c in enumerate(carrier_list)}
            profile = tuple(
                min(cycle_distance(x, c, size) + costs[c] for c in carrier)
                for x in range(size)
            )
            weighted_carriers.append((carrier, costs, profile))

    pair_checks = 0
    exact_zero_cost_checks = 0
    maximum_slack = 0
    for index, (carrier, costs, profile) in enumerate(weighted_carriers):
        for other, other_costs, other_profile in weighted_carriers[index + 1 :]:
            response = max(abs(a - b) for a, b in zip(profile, other_profile))
            directed = max(
                max(min(cycle_distance(c, d, size) for d in other) for c in carrier),
                max(min(cycle_distance(d, c, size) for c in carrier) for d in other),
            )
            slack = abs(response - directed)
            if slack > 1:
                raise AssertionError((carrier, costs, other, other_costs, response, directed))
            maximum_slack = max(maximum_slack, slack)
            if all(value == 0 for value in costs.values()) and all(
                value == 0 for value in other_costs.values()
            ):
                if response != directed:
                    raise AssertionError((carrier, other, response, directed))
                exact_zero_cost_checks += 1
            pair_checks += 1

    return {
        "weighted_carriers": len(weighted_carriers),
        "pair_checks": pair_checks,
        "zero_cost_isometry_checks": exact_zero_cost_checks,
        "maximum_slack_with_p_one": maximum_slack,
    }


def gf16_multiply(a: int, b: int) -> int:
    """Multiply in F_2[x]/(x^4+x+1)."""

    result = 0
    for _ in range(4):
        if b & 1:
            result ^= a
        b >>= 1
        carry = a & 0b1000
        a = (a << 1) & 0b1111
        if carry:
            a ^= 0b0011
    return result


def multiplication_matrix(a: int) -> int:
    """Pack the four columns of multiplication by a into a 16-bit integer."""

    packed = 0
    for column in range(4):
        packed |= gf16_multiply(a, 1 << column) << (4 * column)
    return packed


def rank_4(matrix: int) -> int:
    columns = [(matrix >> (4 * j)) & 0b1111 for j in range(4)]
    rank = 0
    for pivot in range(4):
        pivot_index = next((j for j in range(rank, 4) if (columns[j] >> pivot) & 1), None)
        if pivot_index is None:
            continue
        columns[rank], columns[pivot_index] = columns[pivot_index], columns[rank]
        for j in range(4):
            if j != rank and ((columns[j] >> pivot) & 1):
                columns[j] ^= columns[rank]
        rank += 1
    return rank


def rank_profile(basis: tuple[int, ...], query: int) -> int:
    best = rank_4(query)
    for mask in range(1 << len(basis)):
        shortcut = 0
        cost = 0
        for j, value in enumerate(basis):
            if (mask >> j) & 1:
                shortcut ^= value
                cost += 1
        best = min(best, cost + rank_4(query ^ shortcut))
    return best


def rank_metric_check() -> dict[str, int]:
    field_maps = tuple(multiplication_matrix(a) for a in range(16))
    nonzero_ranks = {rank_4(field_maps[a]) for a in range(1, 16)}
    if nonzero_ranks != {4}:
        raise AssertionError(nonzero_ranks)

    coefficient_subspaces: dict[frozenset[int], tuple[int, int]] = {}
    for a, b in itertools.combinations(range(1, 16), 2):
        coefficient_subspaces.setdefault(span((a, b)), (a, b))
    if len(coefficient_subspaces) != 35:  # [4 choose 2]_2.
        raise AssertionError(len(coefficient_subspaces))

    records = []
    for coefficients, basis_coefficients in coefficient_subspaces.items():
        basis_maps = tuple(field_maps[a] for a in basis_coefficients)
        carrier = frozenset(field_maps[a] for a in coefficients)
        records.append((coefficients, basis_maps, carrier))

    pair_checks = 0
    minimum_carrier_distance = 4
    minimum_response_witness_gap = 4
    for index, (coefficients, basis, carrier) in enumerate(records):
        for other_coefficients, other_basis, other_carrier in records[index + 1 :]:
            witness_coefficient = min(coefficients - other_coefficients)
            witness = field_maps[witness_coefficient]
            carrier_distance = min(rank_4(witness ^ point) for point in other_carrier)
            gap = rank_profile(other_basis, witness) - rank_profile(basis, witness)
            if carrier_distance != 4 or gap < 2:  # D-k = 4-2.
                raise AssertionError((coefficients, other_coefficients, carrier_distance, gap))
            minimum_carrier_distance = min(minimum_carrier_distance, carrier_distance)
            minimum_response_witness_gap = min(minimum_response_witness_gap, gap)
            pair_checks += 1

    # Check CC.15 over every matrix query for a representative subspace.
    _, representative_basis, representative_carrier = records[0]
    profile_formula_checks = 0
    for query in range(1 << 16):
        profile = rank_profile(representative_basis, query)
        carrier_distance = min(rank_4(query ^ point) for point in representative_carrier)
        if not carrier_distance <= profile <= carrier_distance + 2:
            raise AssertionError((query, carrier_distance, profile))
        profile_formula_checks += 1

    return {
        "field_size": 16,
        "matrix_dimension": 4,
        "nonzero_multiplication_rank": min(nonzero_ranks),
        "two_subspaces": len(records),
        "subspace_pairs": pair_checks,
        "minimum_carrier_hausdorff_witness": minimum_carrier_distance,
        "minimum_response_witness_gap": minimum_response_witness_gap,
        "predicted_gap_D_minus_k": 2,
        "profile_formula_checks": profile_formula_checks,
    }


def main() -> None:
    result = {
        "status": "passed",
        "finite_metric": finite_metric_check(),
        "rank_metric": rank_metric_check(),
    }
    RESULT_PATH.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
