#!/usr/bin/env python3
"""Exact checks for contextual block refinement and defect-cycle drift.

The four-state matrices are a nonlinear lift of two-state max-plus matrices.
The exact lift refines contextually to two blocks.  Increasing one microscopic
self-loop by ``delta`` breaks that quotient, refines all four states apart,
and creates worst-case horizon-``n`` response error exactly ``n*delta``.

All calculations use ``Fraction``; no numerical tolerance is involved.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product
import json


Q = Fraction
NEGATIVE_INFINITY = None


def maxplus_vector_matrix(
    vector: tuple[Fraction | None, ...],
    matrix: tuple[tuple[Fraction, ...], ...],
) -> tuple[Fraction | None, ...]:
    result: list[Fraction | None] = []
    for target in range(len(matrix[0])):
        candidates = [
            value + matrix[source][target]
            for source, value in enumerate(vector)
            if value is not NEGATIVE_INFINITY
        ]
        result.append(max(candidates) if candidates else NEGATIVE_INFINITY)
    return tuple(result)


def maxplus_matrix_vector(
    matrix: tuple[tuple[Fraction, ...], ...],
    vector: tuple[Fraction, ...],
) -> tuple[Fraction, ...]:
    return tuple(
        max(matrix[source][target] + vector[target] for target in range(len(vector)))
        for source in range(len(matrix))
    )


def terminal_value(vector: tuple[Fraction | None, ...]) -> Fraction:
    finite = [value for value in vector if value is not NEGATIVE_INFINITY]
    return max(finite)


BLOCKS = ((0, 1), (2, 3))
STATE_BLOCK = (0, 0, 1, 1)
QUOTIENT_A = ((Q(0), Q(-2)), (Q(-1), Q(1)))
QUOTIENT_B = ((Q(-1), Q(0)), (Q(2), Q(-2)))


def lift(quotient: tuple[tuple[Fraction, ...], ...]) -> tuple[tuple[Fraction, ...], ...]:
    rows = []
    for source, source_block in enumerate(STATE_BLOCK):
        microscopic_pattern = (Q(0), Q(-1)) if source in (0, 2) else (Q(-2), Q(0))
        row = []
        for target_block in range(2):
            row.extend(
                quotient[source_block][target_block] + offset
                for offset in microscopic_pattern
            )
        rows.append(tuple(row))
    return tuple(rows)


RAW_A = lift(QUOTIENT_A)
RAW_B = lift(QUOTIENT_B)


def aggregate(vector: tuple[Fraction | None, ...]) -> tuple[Fraction | None, ...]:
    return tuple(
        max(
            value
            for state in block
            if (value := vector[state]) is not NEGATIVE_INFINITY
        )
        if any(vector[state] is not NEGATIVE_INFINITY for state in block)
        else NEGATIVE_INFINITY
        for block in BLOCKS
    )


def partition_from_labels(labels: list[int]) -> set[frozenset[int]]:
    blocks: dict[int, set[int]] = {}
    for state, label in enumerate(labels):
        blocks.setdefault(label, set()).add(state)
    return {frozenset(block) for block in blocks.values()}


def strong_block_refinement(
    matrices: tuple[tuple[tuple[Fraction, ...], ...], ...],
    terminal: tuple[Fraction, ...],
) -> tuple[list[int], list[set[frozenset[int]]]]:
    """Coarsest zero-gauge strong block-max partition."""

    terminal_labels: dict[Fraction, int] = {}
    labels = []
    for value in terminal:
        if value not in terminal_labels:
            terminal_labels[value] = len(terminal_labels)
        labels.append(terminal_labels[value])

    history = [partition_from_labels(labels)]
    while True:
        target_blocks: dict[int, list[int]] = {}
        for state, label in enumerate(labels):
            target_blocks.setdefault(label, []).append(state)
        ordered_labels = sorted(target_blocks)

        signatures = []
        for source in range(len(terminal)):
            transition_signature = tuple(
                max(matrix[source][target] for target in target_blocks[label])
                for matrix in matrices
                for label in ordered_labels
            )
            signatures.append((labels[source], transition_signature))

        signature_labels: dict[tuple[object, ...], int] = {}
        refined = []
        for signature in signatures:
            if signature not in signature_labels:
                signature_labels[signature] = len(signature_labels)
            refined.append(signature_labels[signature])
        if refined == labels:
            return labels, history
        labels = refined
        history.append(partition_from_labels(labels))


def suffix_vectors(
    matrices: tuple[tuple[tuple[Fraction, ...], ...], ...],
    maximum_depth: int,
) -> dict[tuple[int, ...], tuple[Fraction, ...]]:
    vectors: dict[tuple[int, ...], tuple[Fraction, ...]] = {
        (): (Q(0),) * len(matrices[0])
    }
    frontier = [()]
    for _ in range(maximum_depth):
        next_frontier = []
        for word in frontier:
            for letter, matrix in enumerate(matrices):
                extended = (letter,) + word
                vectors[extended] = maxplus_matrix_vector(matrix, vectors[word])
                next_frontier.append(extended)
        frontier = next_frontier
    return vectors


def contextual_partition(
    matrices: tuple[tuple[tuple[Fraction, ...], ...], ...],
    maximum_depth: int,
) -> set[frozenset[int]]:
    vectors = suffix_vectors(matrices, maximum_depth)
    words = sorted(vectors, key=lambda word: (len(word), word))
    signatures: dict[tuple[Fraction, ...], set[int]] = {}
    for state in range(len(matrices[0])):
        signature = tuple(vectors[word][state] for word in words)
        signatures.setdefault(signature, set()).add(state)
    return {frozenset(block) for block in signatures.values()}


def block_defect(
    raw: tuple[tuple[Fraction, ...], ...],
    quotient: tuple[tuple[Fraction, ...], ...],
) -> Fraction:
    return max(
        abs(
            max(raw[source][target] for target in BLOCKS[target_block])
            - quotient[STATE_BLOCK[source]][target_block]
        )
        for source in range(4)
        for target_block in range(2)
    )


def run_forward(
    vector: tuple[Fraction | None, ...],
    word: tuple[int, ...],
    matrices: tuple[tuple[tuple[Fraction, ...], ...], ...],
) -> tuple[Fraction | None, ...]:
    for letter in word:
        vector = maxplus_vector_matrix(vector, matrices[letter])
    return vector


def verify_exact_four_to_two() -> dict[str, object]:
    terminal = (Q(0),) * 4
    labels, history = strong_block_refinement((RAW_A, RAW_B), terminal)
    expected_partition = {frozenset((0, 1)), frozenset((2, 3))}
    assert partition_from_labels(labels) == expected_partition
    assert history == [
        {frozenset((0, 1, 2, 3))},
        expected_partition,
    ]
    assert contextual_partition((RAW_A, RAW_B), 1) == expected_partition
    assert contextual_partition((RAW_A, RAW_B), 4) == expected_partition
    assert block_defect(RAW_A, QUOTIENT_A) == 0
    assert block_defect(RAW_B, QUOTIENT_B) == 0

    # Exact aggregation commutes with both letters on an exhaustive box.
    update_checks = 0
    for entries in product(range(-3, 4), repeat=4):
        raw_vector = tuple(map(Q, entries))
        quotient_vector = aggregate(raw_vector)
        for raw_matrix, quotient_matrix in (
            (RAW_A, QUOTIENT_A),
            (RAW_B, QUOTIENT_B),
        ):
            raw_updated = maxplus_vector_matrix(raw_vector, raw_matrix)
            quotient_updated = maxplus_vector_matrix(
                quotient_vector, quotient_matrix
            )
            assert aggregate(raw_updated) == quotient_updated
            update_checks += 1

    return {
        "coarsest_strong_blocks": [[0, 1], [2, 3]],
        "contextual_depth_one_blocks": [[0, 1], [2, 3]],
        "exhaustive_aggregation_updates": update_checks,
    }


def verify_nearby_failed_quotient() -> dict[str, object]:
    delta = Q(1, 10)
    perturbed_rows = [list(row) for row in RAW_A]
    perturbed_rows[2][2] += delta
    perturbed_a = tuple(tuple(row) for row in perturbed_rows)
    raw_matrices = (perturbed_a, RAW_B)
    quotient_matrices = (QUOTIENT_A, QUOTIENT_B)

    terminal = (Q(0),) * 4
    labels, history = strong_block_refinement(raw_matrices, terminal)
    singleton_partition = {frozenset((state,)) for state in range(4)}
    assert partition_from_labels(labels) == singleton_partition
    assert history == [
        {frozenset((0, 1, 2, 3))},
        {frozenset((0, 1)), frozenset((2,)), frozenset((3,))},
        singleton_partition,
    ]
    assert contextual_partition(raw_matrices, 1) == history[1]
    assert contextual_partition(raw_matrices, 2) == singleton_partition

    assert block_defect(perturbed_a, QUOTIENT_A) == delta
    assert block_defect(RAW_B, QUOTIENT_B) == 0

    # Exhaustively check the word-specific upper bound (# perturbed letters)*delta
    # on a finite box, and record the largest observed error by word length.
    maximum_error_by_length: dict[int, Fraction] = {}
    upper_bound_checks = 0
    words = [word for length in range(6) for word in product((0, 1), repeat=length)]
    for entries in product((-2, 0, 2), repeat=4):
        raw_start = tuple(map(Q, entries))
        quotient_start = aggregate(raw_start)
        for word in words:
            raw_value = terminal_value(run_forward(raw_start, word, raw_matrices))
            quotient_value = terminal_value(
                run_forward(quotient_start, word, quotient_matrices)
            )
            error = abs(raw_value - quotient_value)
            assert error <= word.count(0) * delta
            maximum_error_by_length[len(word)] = max(
                maximum_error_by_length.get(len(word), Q(0)), error
            )
            upper_bound_checks += 1

    # The microscopic cycle 2 --A--> 2 and quotient cycle 1 --A--> 1
    # are both maximizing.  Their edge weights differ by delta, so A^n
    # attains the n*delta upper bound exactly.
    raw_seed = (NEGATIVE_INFINITY, NEGATIVE_INFINITY, Q(0), NEGATIVE_INFINITY)
    quotient_seed = aggregate(raw_seed)
    exact_cycle_checks = 0
    for length in range(1, 41):
        word = (0,) * length
        raw_value = terminal_value(run_forward(raw_seed, word, raw_matrices))
        quotient_value = terminal_value(
            run_forward(quotient_seed, word, quotient_matrices)
        )
        assert raw_value == length * (Q(1) + delta)
        assert quotient_value == Q(length)
        assert raw_value - quotient_value == length * delta
        exact_cycle_checks += 1

    return {
        "perturbed_entry": "T_A[2,2] += 1/10",
        "strong_refinement_block_counts": [len(stage) for stage in history],
        "contextual_singletons_by_depth": 2,
        "uniform_one_step_block_defect": str(delta),
        "exhaustive_upper_bound_checks": upper_bound_checks,
        "maximum_error_by_length_through_five": {
            str(length): str(error)
            for length, error in sorted(maximum_error_by_length.items())
        },
        "exact_linear_cycle_checks": exact_cycle_checks,
        "length_40_cycle_error": str(40 * delta),
    }


def main() -> None:
    print(
        json.dumps(
            {
                "exact_four_to_two": verify_exact_four_to_two(),
                "nearby_failed_quotient": verify_nearby_failed_quotient(),
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
