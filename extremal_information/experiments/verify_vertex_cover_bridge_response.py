#!/usr/bin/env python3
"""Exact finite checks for the vertex-cover bilinear-bridge benchmark."""

from __future__ import annotations

import itertools
import json
import random
from fractions import Fraction


def spins(size):
    return tuple(itertools.product((-1, 1), repeat=size))


def bilinear(left, matrix, right):
    return sum(
        left[i] * matrix[i][j] * right[j]
        for i in range(len(left))
        for j in range(len(right))
    )


def direct_message(landscape, matrix):
    right_states = spins(len(matrix[0]))
    return {
        right: max(
            value + bilinear(left, matrix, right)
            for left, value in landscape.items()
        )
        for right in right_states
    }


def cover_table(landscape, matrix, left_cover_size, right_cover_size):
    left_size = len(matrix)
    table = {}
    for covered_left in spins(left_cover_size):
        for covered_right in spins(right_cover_size):
            covered_term = sum(
                covered_left[i] * matrix[i][j] * covered_right[j]
                for i in range(left_cover_size)
                for j in range(right_cover_size)
            )
            hidden_best = max(
                landscape[covered_left + hidden_left]
                + sum(
                    hidden_left[i - left_cover_size]
                    * matrix[i][j]
                    * covered_right[j]
                    for i in range(left_cover_size, left_size)
                    for j in range(right_cover_size)
                )
                for hidden_left in spins(left_size - left_cover_size)
            )
            table[covered_left, covered_right] = covered_term + hidden_best
    return table


def envelope(
    table, matrix, left_cover_size, right_cover_size, explicit_right_size=None
):
    right_size = (
        explicit_right_size if explicit_right_size is not None else len(matrix[0])
    )
    output = {}
    for covered_right in spins(right_cover_size):
        for free_right in spins(right_size - right_cover_size):
            output[covered_right + free_right] = max(
                table[covered_left, covered_right]
                + sum(
                    covered_left[i]
                    * matrix[i][right_cover_size + j]
                    * free_right[j]
                    for i in range(left_cover_size)
                    for j in range(right_size - right_cover_size)
                )
                for covered_left in spins(left_cover_size)
            )
    return output


def future_response(message, future):
    return max(message[state] + future[state] for state in message)


def verify_cover_factorization(seed=65003):
    rng = random.Random(seed)
    checks = 0
    for left_cover_size in range(3):
        for right_cover_size in range(3):
            for hidden_left_size in range(3):
                for free_right_size in range(3):
                    left_size = left_cover_size + hidden_left_size
                    right_size = right_cover_size + free_right_size
                    if not left_size or not right_size:
                        continue
                    for _ in range(12):
                        matrix = [
                            [
                                0
                                if i >= left_cover_size and j >= right_cover_size
                                else rng.randrange(-3, 4)
                                for j in range(right_size)
                            ]
                            for i in range(left_size)
                        ]
                        landscape = {
                            state: rng.randrange(-8, 9) for state in spins(left_size)
                        }
                        direct = direct_message(landscape, matrix)
                        table = cover_table(
                            landscape, matrix, left_cover_size, right_cover_size
                        )
                        factored = envelope(
                            table, matrix, left_cover_size, right_cover_size
                        )
                        assert direct == factored
                        assert len(table) == 2 ** (
                            left_cover_size + right_cover_size
                        )

                        future = {
                            state: rng.randrange(-9, 10) for state in direct
                        }
                        assert future_response(direct, future) == future_response(
                            factored, future
                        )

                        shift = rng.randrange(-5, 6)
                        shifted_table = {
                            key: value + shift for key, value in table.items()
                        }
                        shifted_message = envelope(
                            shifted_table, matrix, left_cover_size, right_cover_size
                        )
                        assert all(
                            shifted_message[state] == direct[state] + shift
                            for state in direct
                        )
                        checks += len(direct) + 3
    return checks


def pin_future(first, second, target):
    penalty = 1 + max(
        [0]
        + [first[state] - first[target] for state in first]
        + [second[state] - second[target] for state in second]
    )
    return {state: 0 if state == target else -penalty for state in first}


def verify_contextual_pinning(seed=12011):
    rng = random.Random(seed)
    checks = 0
    for size in range(1, 6):
        states = spins(size)
        for _ in range(80):
            first = {state: rng.randrange(-12, 13) for state in states}
            second = {state: rng.randrange(-12, 13) for state in states}
            differences = {state: first[state] - second[state] for state in states}
            for target in states:
                future = pin_future(first, second, target)
                observed = future_response(first, future) - future_response(
                    second, future
                )
                assert observed == differences[target]
                checks += 1
            absolute = max(abs(value) for value in differences.values())
            for _ in range(3):
                future = {state: rng.randrange(-20, 21) for state in states}
                observed = abs(
                    future_response(first, future) - future_response(second, future)
                )
                assert observed <= absolute
                checks += 1
            projective = Fraction(
                max(differences.values()) - min(differences.values()), 2
            )
            midpoint = Fraction(
                max(differences.values()) + min(differences.values()), 2
            )
            assert max(
                abs(Fraction(value) - midpoint) for value in differences.values()
            ) == projective
            checks += 2
    return checks


def matching_matrix(size, weight):
    return [
        [weight if i == j else 0 for j in range(size)] for i in range(size)
    ]


def minimum_matching_cover_size(size):
    # Vertices 0..size-1 are left and size..2*size-1 are right.
    best = 2 * size
    for mask in range(1 << (2 * size)):
        population = bin(mask).count("1")
        if population >= best:
            continue
        if all(
            (mask & (1 << edge)) or (mask & (1 << (size + edge)))
            for edge in range(size)
        ):
            best = population
    return best


def verify_matching_selector(seed=314159):
    rng = random.Random(seed)
    checks = 0
    weight = 5
    spread = 7
    assert 2 * weight > spread
    for size in range(1, 7):
        matrix = matching_matrix(size, weight)
        states = spins(size)
        assert minimum_matching_cover_size(size) == size
        checks += 1
        for _ in range(40):
            landscape = {state: rng.randrange(-spread, 1) for state in states}
            landscape[rng.choice(states)] = 0
            message = direct_message(landscape, matrix)
            assert all(
                message[state] == size * weight + landscape[state]
                for state in states
            )
            checks += len(states)
    return checks


def verify_lattice_count():
    checks = 0
    for cover_size in range(1, 4):
        coordinate_count = 2**cover_size
        for levels in range(1, 4):
            observed = sum(
                max(table) == 0
                for table in itertools.product(
                    range(-levels, 1), repeat=coordinate_count
                )
            )
            predicted = (levels + 1) ** coordinate_count - levels**coordinate_count
            assert observed == predicted
            checks += observed
    return checks


def verify_envelope_nonexpansiveness(seed=424242):
    rng = random.Random(seed)
    checks = 0
    error = Fraction(1, 5)
    for left_cover_size in range(3):
        for right_cover_size in range(3):
            if not left_cover_size and not right_cover_size:
                continue
            right_size = right_cover_size + 2
            matrix = [
                [rng.randrange(-4, 5) for _ in range(right_size)]
                for _ in range(left_cover_size)
            ]
            keys = [
                (left, right)
                for left in spins(left_cover_size)
                for right in spins(right_cover_size)
            ]
            for _ in range(80):
                table = {
                    key: Fraction(rng.randrange(-30, 31), 3) for key in keys
                }
                perturbed = {
                    key: value + rng.choice((-error, Fraction(0), error))
                    for key, value in table.items()
                }
                first = envelope(
                    table,
                    matrix,
                    left_cover_size,
                    right_cover_size,
                    explicit_right_size=right_size,
                )
                second = envelope(
                    perturbed,
                    matrix,
                    left_cover_size,
                    right_cover_size,
                    explicit_right_size=right_size,
                )
                assert max(abs(first[s] - second[s]) for s in first) <= error
                future = {
                    state: Fraction(rng.randrange(-20, 21), 7) for state in first
                }
                assert abs(
                    future_response(first, future)
                    - future_response(second, future)
                ) <= error
                checks += len(first) + 1
    return checks


def main():
    print(
        json.dumps(
            {
                "contextual_pinning_checks": verify_contextual_pinning(),
                "cover_factorization_checks": verify_cover_factorization(),
                "envelope_nonexpansiveness_checks": verify_envelope_nonexpansiveness(),
                "lattice_response_classes_counted": verify_lattice_count(),
                "matching_selector_checks": verify_matching_selector(),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
