#!/usr/bin/env python3
"""Finite audit of directed-margin response amplification examples."""

from __future__ import annotations

import json
from itertools import product
from pathlib import Path


def wt(x: int) -> int:
    return bin(x).count("1")


def line_directed(v: int, w: int, norm) -> int:
    return min(norm(v), norm(v ^ w))


def line_profile(query: int, vector: int, norm) -> int:
    return min(norm(query), 2 + norm(query ^ vector))


def simplex_words() -> tuple[int, ...]:
    words = []
    nonzero_points = tuple(range(1, 8))
    for a in range(1, 8):
        word = 0
        for coordinate, x in enumerate(nonzero_points):
            dot = wt(a & x) % 2
            word |= dot << coordinate
        words.append(word)
    return tuple(words)


def gf8_mul(a: int, b: int) -> int:
    out = 0
    left, right = a, b
    while right:
        if right & 1:
            out ^= left
        right >>= 1
        left <<= 1
        if left & 0b1000:
            left ^= 0b1011  # x^3+x+1
    return out & 0b111


def multiplication_matrix(a: int) -> int:
    # Nine-bit row/column-agnostic encoding; columns are a*1,a*x,a*x^2.
    matrix = 0
    for column in range(3):
        image = gf8_mul(a, 1 << column)
        for row in range(3):
            matrix |= ((image >> row) & 1) << (3 * row + column)
    return matrix


def binary_matrix_rank(encoded: int) -> int:
    rows = [(encoded >> (3 * row)) & 0b111 for row in range(3)]
    rank = 0
    for column in range(3):
        pivot = next((i for i in range(rank, 3) if (rows[i] >> column) & 1), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        for i in range(3):
            if i != rank and ((rows[i] >> column) & 1):
                rows[i] ^= rows[rank]
        rank += 1
    return rank


def product_profile(query, word, vectors, norm) -> int:
    return sum(line_profile(x, vectors[a], norm) for x, a in zip(query, word))


def audit_alphabet(vectors, universe, norm, expected_gap: int) -> dict[str, object]:
    directed = [
        [0 if i == j else line_directed(v, w, norm) for j, w in enumerate(vectors)]
        for i, v in enumerate(vectors)
    ]
    assert all(
        directed[i][j] == (0 if i == j else expected_gap)
        for i in range(7)
        for j in range(7)
    )

    local_response = []
    difference_ranges = {}
    for i in range(7):
        for j in range(7):
            differences = tuple(
                line_profile(x, vectors[i], norm) - line_profile(x, vectors[j], norm)
                for x in universe
            )
            difference_ranges[i, j] = (min(differences), max(differences))
            if i < j:
                distance = max(abs(min(differences)), abs(max(differences)))
                assert distance >= expected_gap - 2
                local_response.append(distance)
    directed_response_table = [
        [difference_ranges[i, j][1] for j in range(7)] for i in range(7)
    ]
    assert all(
        directed_response_table[i][j] == (0 if i == j else expected_gap - 2)
        for i in range(7)
        for j in range(7)
    )

    # Check exact additivity of both directed carrier distances and the
    # response lower bound on every pair of length-two alphabet words.
    words = tuple(product(range(7), repeat=2))
    checked = 0
    for i, left in enumerate(words):
        for right in words[i + 1 :]:
            symbol_distance = sum(a != b for a, b in zip(left, right))
            directed_sum = sum(
                line_directed(vectors[a], vectors[b], norm)
                for a, b in zip(left, right)
            )
            assert directed_sum == expected_gap * symbol_distance
            minimum = sum(difference_ranges[a, b][0] for a, b in zip(left, right))
            maximum = sum(difference_ranges[a, b][1] for a, b in zip(left, right))
            response_distance = max(abs(minimum), abs(maximum))
            assert response_distance >= expected_gap * symbol_distance - 4
            checked += 1
    return {
        "directed_gap": expected_gap,
        "directed_response_table": directed_response_table,
        "local_response_distances": sorted(set(local_response)),
        "length_two_word_pairs_checked": checked,
    }


def main() -> None:
    hamming_vectors = simplex_words()
    assert len(set(hamming_vectors)) == 7
    assert all(wt(v) == 4 for v in hamming_vectors)
    hamming = audit_alphabet(hamming_vectors, tuple(range(1 << 7)), wt, 4)

    rank_vectors = tuple(multiplication_matrix(a) for a in range(1, 8))
    assert len(set(rank_vectors)) == 7
    assert all(binary_matrix_rank(v) == 3 for v in rank_vectors)
    rank = audit_alphabet(
        rank_vectors,
        tuple(range(1 << 9)),
        binary_matrix_rank,
        3,
    )

    result = {
        "status": "passed",
        "binary_simplex_lines": hamming,
        "rank_multiplication_lines": rank,
    }
    output = Path(__file__).with_name(
        "phase3_finite_alphabet_response_amplification.json"
    )
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
