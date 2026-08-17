#!/usr/bin/env python3
"""Exact finite checks for equal word spectra but full path-lift memory."""

from fractions import Fraction as Q
from itertools import product


def maxplus_product(left, right):
    size = len(left)
    return tuple(
        tuple(max(left[i][k] + right[k][j] for k in range(size))
              for j in range(size))
        for i in range(size)
    )


def matrix_for_letter(size, letter, gap):
    return tuple(
        tuple(Q(0) if source == letter else -gap for _ in range(size))
        for source in range(size)
    )


def debruijn_matrix(states, letter, gap):
    return tuple(
        tuple(
            Q(0)
            if source[0] == letter and target[:-1] == source[1:]
            else -gap
            for target in states
        )
        for source in states
    )


def main():
    word_checks = 0
    row_gap_checks = 0
    for size in range(2, 7):
        gap = Q(3, 2)
        matrices = tuple(matrix_for_letter(size, letter, gap)
                         for letter in range(size))
        for depth in range(1, 6):
            for word in product(range(size), repeat=depth):
                product_matrix = matrices[word[0]]
                for letter in word[1:]:
                    product_matrix = maxplus_product(
                        product_matrix, matrices[letter]
                    )
                # Following successive letter names, cyclically, is a
                # zero-weight diagonal path in the word product.
                assert product_matrix[word[0]][word[0]] == 0
                assert max(max(row) for row in product_matrix) == 0
                word_checks += 1

        # Any merged pair has a block-row maximum gap C under one of its two
        # letter names; checking every possible target block is unnecessary
        # here because rows are constant in the target coordinate.
        for left in range(size):
            for right in range(left + 1, size):
                for letter in (left, right):
                    difference = abs(
                        matrices[letter][left][0]
                        - matrices[letter][right][0]
                    )
                    assert difference == gap
                    row_gap_checks += 1

    print(f"zero word-spectrum checks: {word_checks}")
    print(f"full path-lift row-gap checks: {row_gap_checks}")
    print("scalar spectra are exact while every sub-gap path lift is injective")

    binary_checks = 0
    for memory in range(1, 7):
        states = tuple(product((0, 1), repeat=memory))

        def successors(state, letter):
            if state[0] != letter:
                return ()
            return (state[1:] + (0,), state[1:] + (1,))

        def follow_relation(state, word):
            reached = {state}
            for letter in word:
                reached = {
                    target
                    for source in reached
                    for target in successors(source, letter)
                }
            return reached

        # Every binary word has a periodic-window closed lift.
        for depth in range(1, 7):
            for word in product((0, 1), repeat=depth):
                periodic = tuple(word[index % depth]
                                 for index in range(memory))
                assert periodic in follow_relation(periodic, word)
                binary_checks += 1

        # The length-m word naming a state has a zero row only there.
        for state in states:
            for other in states:
                reached = follow_relation(other, state)
                assert bool(reached) == (other == state)
                binary_checks += 1

    print(f"binary de-Bruijn periodic/unique-row checks: {binary_checks}")
    print("a fixed two-letter exact scalar spectrum can require 2^m path-lift states")

    scrambling_checks = 0
    for memory in range(1, 5):
        states = tuple(product((0, 1), repeat=memory))
        gap = Q(memory)
        matrices = tuple(
            debruijn_matrix(states, letter, gap) for letter in (0, 1)
        )
        depth = 2 * memory
        for word in product((0, 1), repeat=depth):
            product_matrix = matrices[word[0]]
            for letter in word[1:]:
                product_matrix = maxplus_product(
                    product_matrix, matrices[letter]
                )
            critical = tuple(word[:memory])
            for state, row in zip(states, product_matrix):
                expected = Q(0) if state == critical else -gap
                assert row == (expected,) * len(states)
                scrambling_checks += 1
            assert product_matrix[states.index(critical)][states.index(critical)] == 0
    print(f"rank-one blocked de-Bruijn checks: {scrambling_checks}")
    print("maximal wordwise contraction still retains exponential rooted state")


if __name__ == "__main__":
    main()
