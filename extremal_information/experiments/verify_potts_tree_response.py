#!/usr/bin/env python3
"""Exact finite checks for the Potts-tree separator-response benchmark.

All score comparisons use integers or fractions.  Randomized coverage is
deterministic under fixed seeds; it is intended to catch algebraic or indexing
errors, while the adjacent benchmark draft contains the general proofs.
"""

from __future__ import annotations

import itertools
import json
import random
from fractions import Fraction


def potts_transfer(predecessor, strength):
    """Max-sum message through K * 1{x == separator}."""
    size = len(predecessor)
    return tuple(
        max(
            predecessor[i] + (strength if i == boundary else 0)
            for i in range(size)
        )
        for boundary in range(size)
    )


def normalize(message):
    baseline = max(message)
    return baseline, tuple(value - baseline for value in message)


def formula_state(predecessor, strength):
    maximum = max(predecessor)
    baseline = maximum + strength
    residual = tuple(max(-strength, value - maximum) for value in predecessor)
    return baseline, residual


def reconstruct(state):
    baseline, residual = state
    return tuple(baseline + value for value in residual)


def response(message, future):
    return max(value + continuation for value, continuation in zip(message, future))


def oscillation(vector):
    return max(vector) - min(vector)


def projective_distance(first, second):
    difference = tuple(a - b for a, b in zip(first, second))
    return Fraction(oscillation(difference), 2)


def verify_clipped_image() -> int:
    checks = 0
    for size in range(2, 5):
        for strength in range(5):
            values = range(-strength - 2, strength + 3)
            for predecessor in itertools.product(values, repeat=size):
                direct = normalize(potts_transfer(predecessor, strength))
                predicted = formula_state(predecessor, strength)
                assert direct == predicted
                baseline, residual = direct
                assert baseline == max(predecessor) + strength
                assert max(residual) == 0
                assert all(-strength <= value <= 0 for value in residual)
                checks += 1

            carrier = [
                residual
                for residual in itertools.product(
                    range(-strength, 1), repeat=size
                )
                if max(residual) == 0
            ]
            assert len(carrier) == (strength + 1) ** size - strength**size
            for residual in carrier:
                chosen_baseline = 7
                predecessor = tuple(
                    chosen_baseline - strength + value for value in residual
                )
                observed = normalize(potts_transfer(predecessor, strength))
                assert observed == (chosen_baseline, residual)
                checks += 1
    return checks


def verify_binary_clamp() -> int:
    checks = 0
    for strength in range(6):
        for first in range(-8, 9):
            for second in range(-8, 9):
                message = potts_transfer((first, second), strength)
                gap = message[1] - message[0]
                predecessor_gap = second - first
                predicted = min(max(predecessor_gap, -strength), strength)
                assert gap == predicted
                checks += 1
    return checks


def pinning_future(first, second, coordinate):
    advantages = [0]
    for vector in (first, second):
        advantages.extend(
            vector[index] - vector[coordinate]
            for index in range(len(vector))
        )
    penalty = max(advantages) + 1
    return tuple(0 if index == coordinate else -penalty for index in range(len(first)))


def verify_contextual_metrics(seed: int = 41041) -> int:
    rng = random.Random(seed)
    checks = 0
    for size in range(2, 7):
        for _ in range(300):
            first = tuple(rng.randint(-12, 12) for _ in range(size))
            second = tuple(rng.randint(-12, 12) for _ in range(size))
            difference = tuple(a - b for a, b in zip(first, second))
            sup_distance = max(abs(value) for value in difference)

            pinned_differences = []
            for coordinate in range(size):
                future = pinning_future(first, second, coordinate)
                observed = response(first, future) - response(second, future)
                assert observed == difference[coordinate]
                pinned_differences.append(observed)
                checks += 1
            assert max(abs(value) for value in pinned_differences) == sup_distance

            midpoint = Fraction(max(difference) + min(difference), 2)
            calibrated_error = max(
                abs(Fraction(value) - midpoint) for value in difference
            )
            assert calibrated_error == Fraction(oscillation(difference), 2)

            for _ in range(5):
                future = tuple(rng.randint(-20, 20) for _ in range(size))
                observed = abs(response(first, future) - response(second, future))
                assert observed <= sup_distance
                checks += 1
    return checks


def brute_tree_message(parent, edge_strength, unary, root_strength):
    size = len(unary[0])
    node_count = len(parent)
    output = []
    for boundary in range(size):
        best = None
        for assignment in itertools.product(range(size), repeat=node_count):
            score = sum(unary[node][assignment[node]] for node in range(node_count))
            score += sum(
                edge_strength[node]
                for node in range(1, node_count)
                if assignment[node] == assignment[parent[node]]
            )
            if assignment[0] == boundary:
                score += root_strength
            if best is None or score > best:
                best = score
        output.append(best)
    return tuple(output)


def recursive_tree_message(parent, edge_strength, unary, root_strength):
    children = [[] for _ in parent]
    for node in range(1, len(parent)):
        children[parent[node]].append(node)

    def visit(node):
        child_messages = [visit(child) for child in children[node]]
        predecessor = tuple(
            unary[node][state]
            + sum(message[state] for message in child_messages)
            for state in range(len(unary[node]))
        )
        strength = root_strength if node == 0 else edge_strength[node]
        return potts_transfer(predecessor, strength)

    return visit(0)


def quotient_tree_message(parent, edge_strength, unary, root_strength):
    children = [[] for _ in parent]
    for node in range(1, len(parent)):
        children[parent[node]].append(node)

    def visit(node):
        child_states = [visit(child) for child in children[node]]
        predecessor_residual = tuple(
            unary[node][state]
            + sum(residual[state] for _, residual in child_states)
            for state in range(len(unary[node]))
        )
        maximum = max(predecessor_residual)
        strength = root_strength if node == 0 else edge_strength[node]
        baseline = (
            sum(child_baseline for child_baseline, _ in child_states)
            + maximum
            + strength
        )
        residual = tuple(
            max(-strength, value - maximum) for value in predecessor_residual
        )
        return baseline, residual

    return reconstruct(visit(0))


def random_tree(rng, node_count, size, fractional=False):
    parent = [-1] + [rng.randrange(node) for node in range(1, node_count)]
    if fractional:
        quantum = Fraction(1, 2)
        edge_strength = [0] + [
            rng.randrange(7) * quantum for _ in range(1, node_count)
        ]
        root_strength = rng.randrange(7) * quantum
        unary = [
            tuple(Fraction(rng.randrange(-30, 31), 6) for _ in range(size))
            for _ in range(node_count)
        ]
    else:
        edge_strength = [0] + [rng.randrange(5) for _ in range(1, node_count)]
        root_strength = rng.randrange(5)
        unary = [
            tuple(rng.randrange(-5, 6) for _ in range(size))
            for _ in range(node_count)
        ]
    return parent, edge_strength, unary, root_strength


def verify_tree_composition(seed: int = 99173) -> int:
    rng = random.Random(seed)
    checks = 0
    for size in (2, 3):
        for node_count in range(1, 7):
            for _ in range(30):
                tree = random_tree(rng, node_count, size)
                brute = brute_tree_message(*tree)
                recursive = recursive_tree_message(*tree)
                quotient = quotient_tree_message(*tree)
                assert brute == recursive == quotient
                checks += size
    return checks


def verify_nonexpansiveness(seed: int = 8128) -> int:
    rng = random.Random(seed)
    checks = 0
    for size in range(2, 7):
        for _ in range(500):
            strength = rng.randrange(8)
            first = tuple(rng.randrange(-20, 21) for _ in range(size))
            second = tuple(rng.randrange(-20, 21) for _ in range(size))
            before = oscillation(tuple(a - b for a, b in zip(first, second)))
            first_image = potts_transfer(first, strength)
            second_image = potts_transfer(second, strength)
            after = oscillation(
                tuple(a - b for a, b in zip(first_image, second_image))
            )
            assert after <= before
            checks += 1

    first = (0, 0)
    second = (0, 1)
    strength = 5
    before = oscillation(tuple(a - b for a, b in zip(first, second)))
    first_image = potts_transfer(first, strength)
    second_image = potts_transfer(second, strength)
    after = oscillation(tuple(a - b for a, b in zip(first_image, second_image)))
    assert after == before == 1
    return checks + 1


def round_to_grid(value: Fraction, quantum: Fraction) -> Fraction:
    scaled = value / quantum
    lower = scaled.numerator // scaled.denominator
    remainder = scaled - lower
    nearest = lower + (1 if remainder > Fraction(1, 2) else 0)
    return nearest * quantum


def verify_lattice_rounding(seed: int = 271828) -> int:
    rng = random.Random(seed)
    quantum = Fraction(1, 2)
    checks = 0
    for size in (2, 3):
        for node_count in range(1, 7):
            for _ in range(20):
                parent, edge_strength, unary, root_strength = random_tree(
                    rng, node_count, size, fractional=True
                )
                rounded_unary = [
                    tuple(round_to_grid(value, quantum) for value in table)
                    for table in unary
                ]
                assert all(
                    abs(value - rounded) <= quantum / 2
                    for table, rounded_table in zip(unary, rounded_unary)
                    for value, rounded in zip(table, rounded_table)
                )

                original = brute_tree_message(
                    parent, edge_strength, unary, root_strength
                )
                rounded = brute_tree_message(
                    parent, edge_strength, rounded_unary, root_strength
                )
                bound = node_count * quantum / 2
                assert all(abs(a - b) <= bound for a, b in zip(original, rounded))

                baseline, residual = normalize(rounded)
                assert (baseline / quantum).denominator == 1
                assert all((value / quantum).denominator == 1 for value in residual)
                assert max(residual) == 0
                assert all(-root_strength <= value <= 0 for value in residual)

                future = tuple(
                    Fraction(rng.randrange(-40, 41), 5) for _ in range(size)
                )
                response_error = abs(
                    response(original, future) - response(rounded, future)
                )
                assert response_error <= bound

                rounded_quotient = quotient_tree_message(
                    parent, edge_strength, rounded_unary, root_strength
                )
                assert rounded_quotient == rounded
                checks += size + 3
    return checks


def verify_entropy_constructions() -> int:
    checks = 0

    # Fine half-grid points are covered by an integer grid on a chosen face.
    strength = Fraction(4)
    epsilon = Fraction(1)
    fine_values = [Fraction(index, 2) for index in range(-8, 1)]
    for size in range(2, 5):
        for residual in itertools.product(fine_values, repeat=size):
            if max(residual) != 0:
                continue
            zero_coordinate = residual.index(0)
            decoded = tuple(
                Fraction(0)
                if index == zero_coordinate
                else round_to_grid(value, epsilon)
                for index, value in enumerate(residual)
            )
            assert max(decoded) == 0
            assert all(-strength <= value <= 0 for value in decoded)
            assert max(abs(a - b) for a, b in zip(residual, decoded)) <= epsilon
            assert projective_distance(residual, decoded) <= epsilon
            checks += 1

    # A 5-epsilon grid on one face is separated by more than 2 epsilon in
    # projective contextual distance.
    strength = 20
    epsilon = 1
    packing_values = (-20, -15, -10)
    for size in range(2, 5):
        packing = [
            tuple(coordinates) + (0,)
            for coordinates in itertools.product(packing_values, repeat=size - 1)
        ]
        for index, first in enumerate(packing):
            for second in packing[index + 1 :]:
                assert projective_distance(first, second) > 2 * epsilon
                assert max(abs(a - b) for a, b in zip(first, second)) > 2 * epsilon
                checks += 1
        assert all(
            max(point) == 0 and all(-strength <= value <= 0 for value in point)
            for point in packing
        )
    return checks


def main() -> None:
    print(
        json.dumps(
            {
                "binary_clamp_checks": verify_binary_clamp(),
                "clipped_image_and_count_checks": verify_clipped_image(),
                "contextual_metric_checks": verify_contextual_metrics(),
                "entropy_cover_and_packing_checks": verify_entropy_constructions(),
                "lattice_rounding_checks": verify_lattice_rounding(),
                "nonexpansiveness_checks": verify_nonexpansiveness(),
                "tree_composition_boundary_checks": verify_tree_composition(),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
