#!/usr/bin/env python3
"""Finite exact checks for the Ising and weighted-residual benchmarks."""

from __future__ import annotations

import itertools
import json
import math
import random


def direct_chain_update(a: int, d: int, coupling: int, field: int) -> tuple[int, int]:
    old = {-1: a, 1: a + d}
    new = {
        y: field * y + max(old[x] + coupling * x * y for x in (-1, 1))
        for y in (-1, 1)
    }
    return new[-1], new[1] - new[-1]


def formula_chain_update(a: int, d: int, coupling: int, field: int) -> tuple[int, int]:
    transmitted = 0
    if coupling:
        clipped = min(max(d, -2 * abs(coupling)), 2 * abs(coupling))
        transmitted = (1 if coupling > 0 else -1) * clipped
    return (
        a - field + max(d - coupling, coupling),
        2 * field + transmitted,
    )


def verify_chain() -> int:
    checks = 0
    for a in range(-2, 3):
        for d in range(-10, 11):
            for coupling in range(-5, 6):
                for field in range(-5, 6):
                    assert direct_chain_update(a, d, coupling, field) == formula_chain_update(
                        a, d, coupling, field
                    )
                    checks += 1
    return checks


def verify_lookup_profiles(seed: int = 1729) -> int:
    rng = random.Random(seed)
    checks = 0
    for width in range(1, 4):
        states = list(itertools.product((-1, 1), repeat=width))
        for _ in range(10):
            table = {s: rng.randint(-7, 7) for s in states}
            offset = min(table.values())
            weights = {s: table[s] - offset for s in states}
            for boundary in states:
                best = -math.inf
                for hidden in itertools.product((0, 1), repeat=len(states)):
                    energy = offset
                    for bit, target in zip(hidden, states):
                        activation = sum(a * s for a, s in zip(target, boundary))
                        energy += bit * weights[target] * (activation - (width - 1))
                    best = max(best, energy)
                assert best == table[boundary]
                checks += 1
    return checks


def vec_mat(vector, matrix):
    return [
        max(vector[i] + matrix[i][j] for i in range(len(vector)))
        for j in range(len(matrix[0]))
    ]


def mat_vec(matrix, vector):
    return [
        max(matrix[i][j] + vector[j] for j in range(len(vector)))
        for i in range(len(matrix))
    ]


def verify_strict_residual_quotient() -> int:
    # The two prefixes expose different raw forward vectors, but every suffix
    # over this fixed transition alphabet gives the same accepted score.
    transition_a = [[0, -1], [0, -1]]
    transition_b = [[0, -7], [0, -7]]
    alpha = [0, -math.inf]
    beta = [0, 0]
    p_a = vec_mat(alpha, transition_a)
    p_b = vec_mat(alpha, transition_b)
    assert p_a != p_b
    checks = 0
    for length in range(8):
        for word in itertools.product((transition_a, transition_b), repeat=length):
            suffix = beta
            for transition in reversed(word):
                suffix = mat_vec(transition, suffix)
            assert max(x + y for x, y in zip(p_a, suffix)) == max(
                x + y for x, y in zip(p_b, suffix)
            )
            checks += 1
    return checks


def verify_directed_ising_bottleneck() -> int:
    checks = 0
    for first in range(-8, 9):
        for second in range(-8, 9):
            if not first or not second:
                continue
            rows = {
                a: {
                    t: min(-first * a * s - second * s * t for s in (-1, 1))
                    for t in (-1, 1)
                }
                for a in (-1, 1)
            }
            forward = max(rows[1][t] - rows[-1][t] for t in (-1, 1))
            reverse = max(rows[-1][t] - rows[1][t] for t in (-1, 1))
            assert forward == reverse == 2 * min(abs(first), abs(second))
            checks += 1
    return checks


def compose_minplus(left, right):
    return [
        [
            min(left[a][u] + right[u][t] for u in range(len(right)))
            for t in range(len(right))
        ]
        for a in range(len(left))
    ]


def permutation_kernel(permutation: tuple[int, ...], strength: int):
    size = len(permutation)
    return [
        [-strength if t == permutation[a] else 0 for t in range(size)]
        for a in range(size)
    ]


def verify_permutation_potts_bottleneck() -> int:
    checks = 0
    for size in range(2, 6):
        permutations = list(itertools.permutations(range(size)))
        # Deterministic thinning keeps q=5 inexpensive while covering varied
        # cycle types and noncommuting pairs.
        stride = max(1, len(permutations) // 12)
        sample = permutations[::stride]
        for first_perm in sample:
            for second_perm in sample:
                composed_perm = tuple(second_perm[first_perm[a]] for a in range(size))
                for first_strength in range(1, 5):
                    for second_strength in range(1, 5):
                        observed = compose_minplus(
                            permutation_kernel(first_perm, first_strength),
                            permutation_kernel(second_perm, second_strength),
                        )
                        baseline = -max(first_strength, second_strength)
                        bottleneck = min(first_strength, second_strength)
                        predicted = [
                            [
                                baseline
                                - (bottleneck if t == composed_perm[a] else 0)
                                for t in range(size)
                            ]
                            for a in range(size)
                        ]
                        assert observed == predicted
                        for a in range(size):
                            for b in range(size):
                                if a == b:
                                    continue
                                directed = max(
                                    observed[a][t] - observed[b][t]
                                    for t in range(size)
                                )
                                assert directed == bottleneck
                        checks += 1
    return checks


def verify_minplus_exposure(seed: int = 8675309) -> int:
    rng = random.Random(seed)
    checks = 0
    for input_size in range(2, 8):
        for output_size in range(1, 7):
            for _ in range(100):
                kernel = [
                    [rng.randint(-8, 8) for _ in range(input_size)]
                    for _ in range(output_size)
                ]
                first = [rng.randint(-8, 8) for _ in range(input_size)]
                second = [rng.randint(-8, 8) for _ in range(input_size)]
                image_first = [
                    min(row[y] + first[y] for y in range(input_size))
                    for row in kernel
                ]
                image_second = [
                    min(row[y] + second[y] for y in range(input_size))
                    for row in kernel
                ]
                before = max(a - b for a, b in zip(first, second))
                after = max(a - b for a, b in zip(image_first, image_second))
                exposed = {
                    y
                    for x, row in enumerate(kernel)
                    for y in range(input_size)
                    if row[y] + first[y] == image_first[x]
                }
                lower = max(first[y] - second[y] for y in exposed)
                assert lower <= after <= before
                checks += 1
    return checks


def main() -> None:
    print(
        json.dumps(
            {
                "chain_recurrence_checks": verify_chain(),
                "lookup_profile_boundary_checks": verify_lookup_profiles(),
                "fixed_automaton_suffix_checks": verify_strict_residual_quotient(),
                "directed_ising_bottleneck_checks": verify_directed_ising_bottleneck(),
                "permutation_potts_bottleneck_checks": verify_permutation_potts_bottleneck(),
                "minplus_exposure_checks": verify_minplus_exposure(),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
