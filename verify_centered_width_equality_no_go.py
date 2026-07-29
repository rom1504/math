#!/usr/bin/env python3
"""Exact no-go for equality four-lifts of the negative triangle.

Equality in the compressed centered-width range bound forces every
cross block to be -J+2P and every diagonal fibre to have its negative
edges form a matching.  Independent fibre permutations reduce the
three cross permutations to I, I, P.  This script enumerates all
24*10^3 reduced candidates and all 2^11 antipodal spin states.
"""

from __future__ import annotations

from itertools import combinations, permutations, product

import numpy as np


def antipodal_spins(order: int) -> np.ndarray:
    spins = np.ones((1 << (order - 1), order), dtype=np.int64)
    for mask in range(1 << (order - 1)):
        for index in range(1, order):
            spins[mask, index] = 1 - 2 * ((mask >> (index - 1)) & 1)
    return spins


def matchings_four() -> list[tuple[tuple[int, int], ...]]:
    pairs = list(combinations(range(4), 2))
    answer = []
    for mask in range(1 << len(pairs)):
        edges = tuple(
            pairs[index]
            for index in range(len(pairs))
            if (mask >> index) & 1
        )
        if all(sum(vertex in edge for edge in edges) <= 1 for vertex in range(4)):
            answer.append(edges)
    assert len(answer) == 10
    return answer


def cross_block(permutation: tuple[int, ...]) -> np.ndarray:
    block = -np.ones((4, 4), dtype=np.int64)
    for row, column in enumerate(permutation):
        block[row, column] = 1
    assert np.all(block.sum(axis=0) == -2)
    assert np.all(block.sum(axis=1) == -2)
    return block


def diagonal_block(matching: tuple[tuple[int, int], ...]) -> np.ndarray:
    block = np.ones((4, 4), dtype=np.int64) - np.eye(4, dtype=np.int64)
    for first, second in matching:
        block[first, second] = block[second, first] = -1
    return block


def candidate(
    permutation: tuple[int, ...],
    matching_indices: tuple[int, int, int],
    matchings: list[tuple[tuple[int, int], ...]],
) -> np.ndarray:
    matrix = np.zeros((12, 12), dtype=np.int64)
    for fibre, matching_index in enumerate(matching_indices):
        matrix[
            4 * fibre : 4 * fibre + 4,
            4 * fibre : 4 * fibre + 4,
        ] = diagonal_block(matchings[matching_index])

    identity_block = cross_block(tuple(range(4)))
    blocks = ((0, 1, identity_block), (0, 2, identity_block))
    for first, second, block in blocks:
        matrix[4 * first : 4 * first + 4, 4 * second : 4 * second + 4] = block
        matrix[4 * second : 4 * second + 4, 4 * first : 4 * first + 4] = block.T

    block = cross_block(permutation)
    matrix[4:8, 8:12] = block
    matrix[8:12, 4:8] = block.T
    return matrix


def main() -> None:
    spins = antipodal_spins(12)
    matchings = matchings_four()
    best_width = None
    multiplicity = 0
    certificate = None
    candidate_count = 0

    for permutation in permutations(range(4)):
        for matching_indices in product(range(10), repeat=3):
            matrix = candidate(permutation, matching_indices, matchings)
            energies = (
                np.einsum("bi,ij,bj->b", spins, matrix, spins) // 2
            )
            lower = int(energies.min())
            upper = int(energies.max())
            width = (upper - lower) // 2
            candidate_count += 1
            if best_width is None or width < best_width:
                best_width = width
                multiplicity = 1
                certificate = {
                    "permutation": permutation,
                    "matching_indices": matching_indices,
                    "lower_endpoint": lower,
                    "upper_endpoint": upper,
                }
            elif width == best_width:
                multiplicity += 1

    assert candidate_count == 24_000
    assert best_width == 20
    assert certificate == {
        "permutation": (0, 1, 2, 3),
        "matching_indices": (0, 0, 0),
        "lower_endpoint": -14,
        "upper_endpoint": 26,
    }
    print("reduced equality candidates:", candidate_count)
    print("minimum centered width:", best_width)
    print("minimum multiplicity:", multiplicity)
    print("first certificate:", certificate)
    print("equality target 16 is impossible")


if __name__ == "__main__":
    main()
