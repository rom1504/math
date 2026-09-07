#!/usr/bin/env python3
"""Finite checks for the independently derived contextual-state formulas."""

from __future__ import annotations

import itertools
import json
import math
import random


def chain_update_direct(offset: int, d: int, coupling: int, field: int) -> tuple[int, int]:
    old = {1: offset + d, -1: offset}
    new = {
        y: field * y + max(old[x] + coupling * x * y for x in (-1, 1))
        for y in (-1, 1)
    }
    return new[-1], new[1] - new[-1]


def chain_update_formula(offset: int, d: int, coupling: int, field: int) -> tuple[int, int]:
    if coupling == 0:
        transmitted = 0
    else:
        clipped = min(max(d, -2 * abs(coupling)), 2 * abs(coupling))
        transmitted = (1 if coupling > 0 else -1) * clipped
    new_offset = offset - field + max(d - coupling, coupling)
    return new_offset, 2 * field + transmitted


def verify_chain_gap() -> int:
    checks = 0
    for offset in range(-3, 4):
        for d in range(-12, 13):
            for coupling in range(-6, 7):
                for field in range(-6, 7):
                    assert chain_update_direct(
                        offset, d, coupling, field
                    ) == chain_update_formula(offset, d, coupling, field)
                    checks += 1
    return checks


def spins(width: int):
    return list(itertools.product((-1, 1), repeat=width))


def verify_lookup_gadget(seed: int = 1729) -> int:
    """Brute-force all hidden spins for arbitrary tables through width three."""
    rng = random.Random(seed)
    checks = 0
    for width in range(1, 4):
        boundary_states = spins(width)
        for _ in range(12):
            table = {state: rng.randint(-9, 9) for state in boundary_states}
            offset = min(table.values())
            weights = {state: value - offset for state, value in table.items()}
            for boundary in boundary_states:
                best = -math.inf
                for hidden in itertools.product((0, 1), repeat=len(boundary_states)):
                    energy = offset
                    for bit, target in zip(hidden, boundary_states):
                        activation = sum(a * s for a, s in zip(target, boundary))
                        energy += bit * weights[target] * (activation - (width - 1))
                    best = max(best, energy)
                assert best == table[boundary]
                checks += 1
    return checks


def verify_selector_isometry(seed: int = 2718) -> int:
    rng = random.Random(seed)
    checks = 0
    for dimension in range(2, 17):
        for _ in range(50):
            first = [rng.randint(-20, 20) for _ in range(dimension)]
            second = [rng.randint(-20, 20) for _ in range(dimension)]
            penalty = 1 + max(max(first) - min(first), max(second) - min(second))
            observed = []
            for selected in range(dimension):
                context = [0 if j == selected else -penalty for j in range(dimension)]
                left = max(x + z for x, z in zip(first, context))
                right = max(x + z for x, z in zip(second, context))
                observed.append(abs(left - right))
            assert max(observed) == max(abs(x - y) for x, y in zip(first, second))
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


def mat_mat(left, right):
    return [
        [
            max(left[i][k] + right[k][j] for k in range(len(right)))
            for j in range(len(right[0]))
        ]
        for i in range(len(left))
    ]


def verify_weighted_automaton(seed: int = 31415) -> dict[str, int]:
    rng = random.Random(seed)
    associative = 0
    nonexpansive = 0
    prefix_suffix = 0
    for size in range(1, 6):
        for _ in range(100):
            matrices = [
                [[rng.randint(-8, 8) for _ in range(size)] for _ in range(size)]
                for _ in range(3)
            ]
            left = mat_mat(mat_mat(matrices[0], matrices[1]), matrices[2])
            right = mat_mat(matrices[0], mat_mat(matrices[1], matrices[2]))
            assert left == right
            associative += 1

            p = [rng.randint(-8, 8) for _ in range(size)]
            perturbation = [rng.randint(-3, 3) for _ in range(size)]
            p2 = [x + z for x, z in zip(p, perturbation)]
            before = max(abs(x - y) for x, y in zip(p, p2))
            after_vectors = vec_mat(p, matrices[0]), vec_mat(p2, matrices[0])
            after = max(abs(x - y) for x, y in zip(*after_vectors))
            assert after <= before
            nonexpansive += 1

            alpha = [rng.randint(-5, 5) for _ in range(size)]
            beta = [rng.randint(-5, 5) for _ in range(size)]
            prefix = vec_mat(alpha, matrices[0])
            suffix = mat_vec(mat_mat(matrices[1], matrices[2]), beta)
            factored = max(x + y for x, y in zip(prefix, suffix))
            whole = max(
                x + y
                for x, y in zip(
                    vec_mat(alpha, mat_mat(matrices[0], mat_mat(matrices[1], matrices[2]))),
                    beta,
                )
            )
            assert factored == whole
            prefix_suffix += 1

    # Two actual prefixes can have distinct forward vectors but identical
    # responses to every suffix in a fixed automaton.
    transition_a = [[0, -1], [0, -1]]
    transition_b = [[0, -7], [0, -7]]
    alpha = [0, -math.inf]
    beta = [0, 0]
    p_a = vec_mat(alpha, transition_a)
    p_b = vec_mat(alpha, transition_b)
    assert p_a != p_b
    for word_length in range(7):
        for word in itertools.product((transition_a, transition_b), repeat=word_length):
            h = beta
            for transition in reversed(word):
                h = mat_vec(transition, h)
            assert max(x + y for x, y in zip(p_a, h)) == max(
                x + y for x, y in zip(p_b, h)
            )

    return {
        "associativity_checks": associative,
        "nonexpansive_checks": nonexpansive,
        "prefix_suffix_checks": prefix_suffix,
        "fixed_automaton_quotient_suffixes": sum(2**length for length in range(7)),
    }


def main() -> None:
    result = {
        "chain_gap_checks": verify_chain_gap(),
        "lookup_gadget_boundary_checks": verify_lookup_gadget(),
        "selector_isometry_random_pairs": verify_selector_isometry(),
    }
    result.update(verify_weighted_automaton())
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
