#!/usr/bin/env python3
"""Finite checks for the scale-rank response sandwich and MRD host."""

from __future__ import annotations

import itertools
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
RESULT_PATH = HERE / "phase3_scale_rank_response_sandwich_results.json"


def binary_span(basis: tuple[int, ...]) -> frozenset[int]:
    values = {0}
    for vector in basis:
        values |= {value ^ vector for value in tuple(values)}
    return frozenset(values)


def all_binary_subspaces(dimension: int) -> tuple[frozenset[int], ...]:
    subspaces = {frozenset({0})}
    for size in range(1, dimension + 1):
        for basis in itertools.combinations(range(1, 1 << dimension), size):
            span = binary_span(basis)
            if len(span) == 1 << size:
                subspaces.add(span)
    return tuple(subspaces)


def two_scale_rank_check() -> dict[str, object]:
    dimension = 3
    quotient_rank = 1
    large_scale = 4
    subspaces = all_binary_subspaces(dimension)

    def weight(vector: int) -> int:
        if vector == 0:
            return 0
        return large_scale + 1 if vector & 1 else 1

    observed: dict[int, int] = {}
    for threshold in (0, 1, 4, 5):
        best = 0
        for subspace in subspaces:
            nonzero = [weight(vector) for vector in subspace if vector]
            if not nonzero or min(nonzero) > threshold:
                best = max(best, (len(subspace)).bit_length() - 1)
        observed[threshold] = best
    expected = {0: 3, 1: 1, 4: 1, 5: 0}
    if observed != expected:
        raise AssertionError((observed, expected))
    return {
        "ambient_dimension": dimension,
        "quotient_dimension": quotient_rank,
        "scale_ranks": {str(key): value for key, value in observed.items()},
        "subspaces_checked": len(subspaces),
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


def gf16_square(a: int) -> int:
    return gf16_multiply(a, a)


def linearized_matrix(a0: int, a1: int) -> int:
    """Pack x -> a0*x + a1*x^2 as four binary columns."""

    packed = 0
    for column in range(4):
        x = 1 << column
        image = gf16_multiply(a0, x) ^ gf16_multiply(a1, gf16_square(x))
        packed |= image << (4 * column)
    return packed


def rank_four(matrix: int) -> int:
    columns = [(matrix >> (4 * column)) & 0b1111 for column in range(4)]
    rank = 0
    for pivot in range(4):
        pivot_index = next((index for index in range(rank, 4) if (columns[index] >> pivot) & 1), None)
        if pivot_index is None:
            continue
        columns[rank], columns[pivot_index] = columns[pivot_index], columns[rank]
        for index in range(4):
            if index != rank and ((columns[index] >> pivot) & 1):
                columns[index] ^= columns[rank]
        rank += 1
    return rank


def gabidulin_host_check() -> dict[str, int]:
    host = {
        linearized_matrix(a0, a1)
        for a0 in range(16)
        for a1 in range(16)
    }
    if len(host) != 256:
        raise AssertionError(len(host))
    ranks = {matrix: rank_four(matrix) for matrix in host}
    minimum_rank = min(rank for matrix, rank in ranks.items() if matrix)
    if minimum_rank != 3:
        raise AssertionError(minimum_rank)

    nonzero = tuple(matrix for matrix in host if matrix)
    minimum_profile_gap = 4
    pair_checks = 0
    for index, left in enumerate(nonzero):
        left_at_left = min(rank_four(left), 2)
        for right in nonzero[index + 1 :]:
            right_at_left = min(rank_four(left), 2 + rank_four(left ^ right))
            gap = right_at_left - left_at_left
            if gap < 1:  # minimum rank three minus presentation radius two.
                raise AssertionError((left, right, left_at_left, right_at_left))
            minimum_profile_gap = min(minimum_profile_gap, gap)
            pair_checks += 1

    return {
        "matrix_dimension": 4,
        "gabidulin_parameter_r": 2,
        "host_dimension": 8,
        "host_size": len(host),
        "minimum_rank_distance": minimum_rank,
        "one_subspace_profiles": len(nonzero),
        "profile_pairs_checked": pair_checks,
        "minimum_witness_gap": minimum_profile_gap,
        "proved_gap": 1,
    }


def main() -> None:
    result = {
        "status": "passed",
        "two_scale_metric": two_scale_rank_check(),
        "rank_metric_mrd": gabidulin_host_check(),
    }
    RESULT_PATH.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
