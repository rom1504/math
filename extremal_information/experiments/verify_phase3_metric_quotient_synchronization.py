#!/usr/bin/env python3
"""Finite checks for metric-quotient synchronization."""

from __future__ import annotations

import itertools
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
RESULT_PATH = HERE / "phase3_metric_quotient_synchronization_results.json"


def cycle_distance(x: int, y: int, size: int) -> int:
    delta = abs(x - y)
    return min(delta, size - delta)


def product_metric_check() -> dict[str, int]:
    """Exhaust all zero-one presentations on C3 times a two-point fibre."""

    points = tuple(itertools.product(range(3), range(2)))

    def distance(left: tuple[int, int], right: tuple[int, int]) -> int:
        return cycle_distance(left[0], right[0], 3) + abs(left[1] - right[1])

    weighted_carriers = 0
    endpoint_checks = 0
    maximum_error = 0
    for labels in itertools.product((-1, 0, 1), repeat=len(points)):
        carrier = [point for point, label in zip(points, labels) if label >= 0]
        if not carrier:
            continue
        costs = {point: labels[index] for index, point in enumerate(points) if labels[index] >= 0}
        quotient_carrier = {point[0] for point in carrier}
        for query in points:
            response = min(distance(query, point) + costs[point] for point in carrier)
            decoder = min(cycle_distance(query[0], y, 3) for y in quotient_carrier)
            error = response - decoder
            if not 0 <= error <= 2:  # fibre diameter one plus presentation radius one.
                raise AssertionError((carrier, costs, query, response, decoder))
            maximum_error = max(maximum_error, error)
            endpoint_checks += 1
        weighted_carriers += 1
    return {
        "weighted_carriers": weighted_carriers,
        "endpoint_checks": endpoint_checks,
        "maximum_error": maximum_error,
        "proved_bound": 2,
    }


def ternary_vectors(dimension: int) -> tuple[tuple[int, ...], ...]:
    return tuple(itertools.product(range(3), repeat=dimension))


def two_scale_check() -> dict[str, int]:
    """Check the fixed-rank quotient decoder in the two-scale carrier."""

    vectors = ternary_vectors(3)
    nonzero = tuple(vector for vector in vectors if any(vector))
    scale = 5

    def add(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
        return tuple((a + b) % 3 for a, b in zip(left, right))

    def multiply(scalar: int, vector: tuple[int, ...]) -> tuple[int, ...]:
        return tuple((scalar * value) % 3 for value in vector)

    def distance(left: tuple[int, ...], right: tuple[int, ...]) -> int:
        difference = tuple((a - b) % 3 for a, b in zip(left, right))
        return scale * int(difference[0] != 0) + int(any(difference))

    maximum_error = 0
    checks = 0
    quotient_states: set[tuple[int, ...]] = set()
    for column in nonzero:
        carrier = {multiply(scalar, column) for scalar in range(3)}
        projected = tuple(sorted({point[0] for point in carrier}))
        quotient_states.add(projected)
        for query in vectors:
            response = min(
                2 * int(scalar != 0) + distance(query, multiply(scalar, column))
                for scalar in range(3)
            )
            decoder = (scale + 1) * int(query[0] not in projected)
            error = abs(response - decoder)
            if error > 3:  # 2k plus fibre diameter one.
                raise AssertionError((column, query, response, decoder))
            maximum_error = max(maximum_error, error)
            checks += 1
    return {
        "maps": len(nonzero),
        "endpoint_checks": checks,
        "quotient_carrier_states": len(quotient_states),
        "maximum_error": maximum_error,
        "proved_bound": 3,
    }


def rank_two(matrix: int) -> int:
    rows = ((matrix >> 2) & 0b11, matrix & 0b11)
    nonzero = [row for row in rows if row]
    if not nonzero:
        return 0
    if len(nonzero) == 1 or nonzero[0] == nonzero[1]:
        return 1
    return 2


def rank_projection_check() -> dict[str, int]:
    """Exhaust all carriers in M_2(F_2) under first-row projection."""

    points = tuple(range(16))
    maximum_error = 0
    endpoint_checks = 0
    for mask in range(1, 1 << len(points)):
        carrier = tuple(point for point in points if (mask >> point) & 1)
        projected = {(point >> 2) & 0b11 for point in carrier}
        for query in points:
            carrier_distance = min(rank_two(query ^ point) for point in carrier)
            quotient_distance = int(((query >> 2) & 0b11) not in projected)
            error = carrier_distance - quotient_distance
            if not 0 <= error <= 1:  # D-r=1.
                raise AssertionError((carrier, query, carrier_distance, quotient_distance))
            maximum_error = max(maximum_error, error)
            endpoint_checks += 1
    return {
        "carriers": (1 << len(points)) - 1,
        "endpoint_checks": endpoint_checks,
        "maximum_error": maximum_error,
        "proved_bound": 1,
    }


def min_plus_check() -> dict[str, int]:
    """Check nonexpansiveness on a finite exhaustive family of profiles."""

    profiles = tuple(itertools.product(range(3), repeat=4))
    kernels = tuple(
        tuple(tuple((multiplier * x + z + shift) % 4 for z in range(4)) for x in range(4))
        for multiplier in range(4)
        for shift in range(4)
    )
    checks = 0
    for left in profiles:
        for right in profiles:
            error = max(abs(a - b) for a, b in zip(left, right))
            for kernel in kernels:
                transformed_left = tuple(
                    min(left[z] + kernel[x][z] for z in range(4)) for x in range(4)
                )
                transformed_right = tuple(
                    min(right[z] + kernel[x][z] for z in range(4)) for x in range(4)
                )
                transformed_error = max(
                    abs(a - b) for a, b in zip(transformed_left, transformed_right)
                )
                if transformed_error > error:
                    raise AssertionError((left, right, kernel, error, transformed_error))
                checks += 1
    return {"profile_kernel_checks": checks}


def main() -> None:
    result = {
        "status": "passed",
        "product_metric": product_metric_check(),
        "two_scale_metric": two_scale_check(),
        "rank_row_projection": rank_projection_check(),
        "min_plus_nonamplification": min_plus_check(),
    }
    RESULT_PATH.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
