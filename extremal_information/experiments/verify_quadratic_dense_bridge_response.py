#!/usr/bin/env python3
"""Exact checks for the quadratic-child dense-bridge packing theorem."""

from __future__ import annotations

import itertools
import json
import random
from fractions import Fraction


def spins(size):
    return tuple(itertools.product((-1, 1), repeat=size))


def dot(first, second):
    return sum(a * b for a, b in zip(first, second))


def matvec(matrix, vector):
    return tuple(dot(row, vector) for row in matrix)


def sign(vector):
    return tuple(1 if value >= 0 else -1 for value in vector)


def pole_quadratic(state, pole):
    size = len(state)
    overlap = dot(state, pole)
    return (overlap * overlap - size) // 2


def verify_pole_lock(seed=20260817):
    rng = random.Random(seed)
    checks = 0
    for size in range(3, 10):
        states = spins(size)
        coordinate_bound = (size - 1) // 2
        for _ in range(200):
            pole = tuple(rng.choice((-1, 1)) for _ in range(size))
            field = tuple(
                rng.randint(-coordinate_bound, coordinate_bound)
                for _ in range(size)
            )
            assert 2 * max(abs(value) for value in field) < size
            observed = max(
                pole_quadratic(state, pole) + dot(field, state)
                for state in states
            )
            predicted = size * (size - 1) // 2 + abs(dot(field, pole))
            assert observed == predicted
            assert all(
                pole[i] * pole[j] in (-1, 1)
                for i in range(size)
                for j in range(i + 1, size)
            )
            checks += 1
    return checks


def bridge_trial(size, query_count, rng):
    matrix = [
        [rng.choice((-1, 1)) for _ in range(size)] for _ in range(size)
    ]
    queries = [
        tuple(rng.choice((-1, 1)) for _ in range(size))
        for _ in range(query_count)
    ]
    fields = [matvec(matrix, query) for query in queries]
    poles = [sign(field) for field in fields]
    maximum_coordinate = max(abs(value) for field in fields for value in field)
    if 2 * maximum_coordinate >= size:
        return None

    canonical_poles = {
        min(pole, tuple(-value for value in pole)) for pole in poles
    }
    if len(canonical_poles) != query_count:
        return None

    diagonals = [sum(abs(value) for value in field) for field in fields]
    response_matrix = [
        [size * (size - 1) // 2 + abs(dot(pole, field)) for field in fields]
        for pole in poles
    ]
    directed_gaps = [
        response_matrix[child][child] - response_matrix[other][child]
        for child in range(query_count)
        for other in range(query_count)
        if child != other
    ]
    projective_gaps = []
    for first in range(query_count):
        for second in range(first + 1, query_count):
            difference = [
                response_matrix[first][query] - response_matrix[second][query]
                for query in range(query_count)
            ]
            projective_gaps.append(
                Fraction(max(difference) - min(difference), 2)
            )

    return {
        "matrix": matrix,
        "queries": queries,
        "poles": poles,
        "minimum_diagonal": min(diagonals),
        "maximum_coordinate": maximum_coordinate,
        "minimum_directed_gap": min(directed_gaps),
        "minimum_projective_gap": min(projective_gaps),
    }


def verify_finite_certificate(seed=8675309, trials=500):
    rng = random.Random(seed)
    size = 32
    query_count = 12
    best = None
    for _ in range(trials):
        candidate = bridge_trial(size, query_count, rng)
        if candidate is None:
            continue
        if best is None or candidate["minimum_directed_gap"] > best[
            "minimum_directed_gap"
        ]:
            best = candidate

    assert best is not None
    assert best["minimum_directed_gap"] > 0
    assert best["minimum_projective_gap"] > 0
    assert 2 * best["maximum_coordinate"] < size
    assert all(
        entry in (-1, 1) for row in best["matrix"] for entry in row
    )
    assert all(
        pole[i] * pole[j] in (-1, 1)
        for pole in best["poles"]
        for i in range(size)
        for j in range(i + 1, size)
    )

    return {
        "n": size,
        "query_count": query_count,
        "trials": trials,
        "minimum_diagonal": best["minimum_diagonal"],
        "maximum_field_coordinate": best["maximum_coordinate"],
        "minimum_absolute_gap": best["minimum_directed_gap"],
        "minimum_projective_gap": float(best["minimum_projective_gap"]),
        "normalized_absolute_gap_n_3_over_2": (
            best["minimum_directed_gap"] / size**1.5
        ),
        "arithmetic": "exact integer except displayed normalization",
    }


def round_to_grid(value, mesh):
    scaled = value / mesh
    lower = scaled.numerator // scaled.denominator
    remainder = scaled - lower
    nearest = lower + (1 if remainder > Fraction(1, 2) else 0)
    return nearest * mesh


def quadratic_value(state, coefficients):
    return sum(
        value * state[first] * state[second]
        for (first, second), value in coefficients.items()
    )


def bridge_message(coefficients, matrix, states):
    landscape = {
        state: quadratic_value(state, coefficients) for state in states
    }
    return {
        query: max(
            landscape[state] + dot(state, matvec(matrix, query))
            for state in states
        )
        for query in states
    }


def future_response(message, future):
    return max(message[state] + future[state] for state in message)


def verify_coefficient_quantization(seed=271828):
    rng = random.Random(seed)
    mesh = Fraction(1, 4)
    checks = 0
    for size in range(2, 8):
        states = spins(size)
        edges = tuple(itertools.combinations(range(size), 2))
        bound = len(edges) * mesh / 2
        for _ in range(12):
            coefficients = {
                edge: Fraction(rng.randint(-20, 20), 20) for edge in edges
            }
            rounded = {
                edge: round_to_grid(value, mesh)
                for edge, value in coefficients.items()
            }
            assert all(
                abs(coefficients[edge] - rounded[edge]) <= mesh / 2
                for edge in edges
            )
            assert max(
                abs(
                    quadratic_value(state, coefficients)
                    - quadratic_value(state, rounded)
                )
                for state in states
            ) <= bound

            matrix = [
                [rng.choice((-1, 1)) for _ in range(size)]
                for _ in range(size)
            ]
            first = bridge_message(coefficients, matrix, states)
            second = bridge_message(rounded, matrix, states)
            assert max(abs(first[state] - second[state]) for state in states) <= bound

            future = {
                state: Fraction(rng.randint(-30, 30), 7) for state in states
            }
            assert abs(
                future_response(first, future) - future_response(second, future)
            ) <= bound
            checks += len(states) + 3
    return checks


def main():
    print(
        json.dumps(
            {
                "coefficient_quantization_checks": verify_coefficient_quantization(),
                "finite_dense_bridge_certificate": verify_finite_certificate(),
                "pole_lock_exhaustive_checks": verify_pole_lock(),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
