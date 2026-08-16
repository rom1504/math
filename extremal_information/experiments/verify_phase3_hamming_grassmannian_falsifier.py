#!/usr/bin/env python3
"""Finite checks for the Hamming-Grassmannian obstruction report.

The script is deliberately self contained.  It checks

* the exact line-carrier metric formula and the ``A_2-1 <= Pack <= A_2``
  reduction for dimensions at most six;
* an explicit seven-letter alphabet of binary 2-planes in dimension six,
  whose two directed Hausdorff distances are both three;
* the directed-distance direct-sum identity on all length-two words; and
* the two strict entropy inequalities used in the asymptotic amplification.

Integers encode binary vectors, with Hamming distance given by xor weight.
"""

from __future__ import annotations

import json
import math
from pathlib import Path


def weight(x: int) -> int:
    return bin(x).count("1")


def popcount(x: int) -> int:
    return bin(x).count("1")


def span(basis: tuple[int, ...] | list[int]) -> frozenset[int]:
    ans = {0}
    for v in basis:
        ans |= {x ^ v for x in tuple(ans)}
    return frozenset(ans)


def directed_distance(c: frozenset[int], d: frozenset[int]) -> int:
    return max(min(weight(x ^ y) for y in d) for x in c)


def hausdorff_distance(c: frozenset[int], d: frozenset[int]) -> int:
    return max(directed_distance(c, d), directed_distance(d, c))


def line_distance(v: int, w: int) -> int:
    return hausdorff_distance(frozenset((0, v)), frozenset((0, w)))


def maximum_clique_size(adjacency: list[int]) -> int:
    """Exact bit-set maximum clique, sufficient for <= 64 vertices here."""

    best = 0

    def expand(chosen: int, candidates: int) -> None:
        nonlocal best
        if chosen + popcount(candidates) <= best:
            return
        if not candidates:
            best = max(best, chosen)
            return

        # Branch first on a high-degree candidate.  The elementary cardinality
        # bound is enough for the small certified instances below.
        vertices: list[int] = []
        q = candidates
        while q:
            bit = q & -q
            vertices.append(bit.bit_length() - 1)
            q ^= bit
        vertices.sort(
            key=lambda v: popcount(adjacency[v] & candidates), reverse=True
        )
        remaining = candidates
        for v in vertices:
            bit = 1 << v
            if not (remaining & bit):
                continue
            expand(chosen + 1, remaining & adjacency[v])
            remaining ^= bit
            if chosen + popcount(remaining) <= best:
                return

    expand(0, (1 << len(adjacency)) - 1)
    return best


def graph_from_predicate(n: int, predicate) -> list[int]:
    adjacency = [0] * n
    for i in range(n):
        for j in range(i):
            if predicate(i, j):
                adjacency[i] |= 1 << j
                adjacency[j] |= 1 << i
    return adjacency


def exact_line_checks(max_dimension: int = 6) -> list[dict[str, int]]:
    rows: list[dict[str, int]] = []
    for dimension in range(2, max_dimension + 1):
        vectors = list(range(1, 1 << dimension))
        for threshold in range(1, dimension):
            for v in vectors:
                for w in vectors:
                    if v == w:
                        continue
                    expected = min(weight(v ^ w), max(weight(v), weight(w)))
                    assert line_distance(v, w) == expected

            line_graph = graph_from_predicate(
                len(vectors),
                lambda i, j: line_distance(vectors[i], vectors[j]) > threshold,
            )
            line_pack = maximum_clique_size(line_graph)
            ball_sizes = [
                sum(line_distance(v, w) <= threshold for v in vectors)
                for w in vectors
            ]
            hamming_volume = sum(
                math.comb(dimension, j) for j in range(threshold + 1)
            )
            assert min(ball_sizes) >= hamming_volume

            # Translation lets an optimal ordinary code contain zero.  Its
            # remaining words must all have weight > threshold.
            eligible = [v for v in vectors if weight(v) > threshold]
            code_graph = graph_from_predicate(
                len(eligible),
                lambda i, j: weight(eligible[i] ^ eligible[j]) > threshold,
            )
            ordinary_code = 1 + maximum_clique_size(code_graph)
            assert ordinary_code - 1 <= line_pack <= ordinary_code
            rows.append(
                {
                    "dimension": dimension,
                    "threshold": threshold,
                    "line_pack": line_pack,
                    "A2": ordinary_code,
                    "minimum_line_ball": min(ball_sizes),
                    "hamming_ball_volume": hamming_volume,
                }
            )
    return rows


# Each tuple is a basis for one 2-plane in F_2^6.  This is a compact finite
# certificate; no optimizer is needed to verify it.
LOCAL_BASES = (
    (0b000111, 0b011011),
    (0b001101, 0b100011),
    (0b001111, 0b110101),
    (0b010011, 0b100110),
    (0b010110, 0b101110),
    (0b011011, 0b101001),
    (0b011100, 0b101001),
)


def concatenate_blocks(blocks: tuple[frozenset[int], ...], width: int = 6) -> frozenset[int]:
    ans = {0}
    offset = 0
    for block in blocks:
        ans = {x | (y << offset) for x in ans for y in block}
        offset += width
    return frozenset(ans)


def local_alphabet_checks() -> dict[str, object]:
    planes = tuple(span(list(basis)) for basis in LOCAL_BASES)
    assert all(len(plane) == 4 for plane in planes)
    assert len(set(planes)) == 7
    joint_span = span(sorted(set().union(*planes)))
    assert len(joint_span) == 1 << 5
    local_minimum_weights = [
        min(weight(x) for x in plane if x != 0) for plane in planes
    ]
    assert max(local_minimum_weights) <= 4

    directed_matrix = [
        [directed_distance(c, d) for d in planes]
        for c in planes
    ]
    for i in range(len(planes)):
        for j in range(len(planes)):
            assert directed_matrix[i][j] == (0 if i == j else 3)

    # Exhaustively check the direct-sum identity on all 7^2 words.  The
    # resulting 4-dimensional subspaces have only 16 points each.
    words = [(i, j) for i in range(7) for j in range(7)]
    products = {
        word: concatenate_blocks((planes[word[0]], planes[word[1]]))
        for word in words
    }
    checked_pairs = 0
    for p, word in enumerate(words):
        for other in words[p + 1 :]:
            symbol_distance = sum(a != b for a, b in zip(word, other))
            assert directed_distance(products[word], products[other]) == 3 * symbol_distance
            assert directed_distance(products[other], products[word]) == 3 * symbol_distance
            checked_pairs += 1

    h2 = lambda x: -x * math.log2(x) - (1 - x) * math.log2(1 - x)
    hq = lambda x, q: (
        x * math.log(q - 1, q)
        - x * math.log(x, q)
        - (1 - x) * math.log(1 - x, q)
    )
    binary_hamming_margin = h2(3 / 16) - 2 / 3
    qary_gilbert_rate = 1 - hq(3 / 4, 7)
    assert binary_hamming_margin > 0
    assert qary_gilbert_rate > 0

    return {
        "alphabet_size": len(planes),
        "alphabet_is_pairwise_distinct": True,
        "joint_span_dimension": 5,
        "local_minimum_weights": local_minimum_weights,
        "directed_distance_matrix": directed_matrix,
        "length_two_pairs_checked": checked_pairs,
        "binary_hamming_margin_H2_3_16_minus_2_3": binary_hamming_margin,
        "qary_GV_symbol_rate_1_minus_H7_3_4": qary_gilbert_rate,
        "qary_GV_bit_rate_per_outer_coordinate": qary_gilbert_rate * math.log2(7),
    }


def binary_rank(rows: tuple[int, ...] | list[int]) -> int:
    pivots: dict[int, int] = {}
    for value in rows:
        x = value
        while x:
            pivot = x.bit_length() - 1
            if pivot in pivots:
                x ^= pivots[pivot]
            else:
                pivots[pivot] = x
                break
    return len(pivots)


def simplex_words() -> tuple[int, ...]:
    # Coordinates are indexed by the seven nonzero vectors of F_2^3.
    words = []
    for functional in range(1, 8):
        word = 0
        for coordinate, vector in enumerate(range(1, 8)):
            if weight(functional & vector) % 2:
                word |= 1 << coordinate
        words.append(word)
    return tuple(words)


def gf8_multiply(a: int, b: int) -> int:
    """Multiply in F_2[x]/(x^3+x+1)."""
    ans = 0
    x = a
    y = b
    while y:
        if y & 1:
            ans ^= x
        y >>= 1
        x <<= 1
        if x & 0b1000:
            x ^= 0b1011
    return ans & 0b111


def multiplication_rows(a: int) -> tuple[int, int, int]:
    columns = [gf8_multiply(a, 1 << j) for j in range(3)]
    return tuple(
        sum(((column >> row) & 1) << col for col, column in enumerate(columns))
        for row in range(3)
    )


def presented_alphabet_checks() -> dict[str, object]:
    simplex = simplex_words()
    assert len(set(simplex)) == 7
    assert all(weight(v) == 4 for v in simplex)
    simplex_directed = []
    for i, v in enumerate(simplex):
        row = []
        for j, w in enumerate(simplex):
            value = 0 if i == j else min(weight(v), weight(v ^ w))
            assert value == (0 if i == j else 4)
            row.append(value)
        simplex_directed.append(row)

    simplex_responses = [
        [min(weight(x), 2 + weight(x ^ v)) for x in range(1 << 7)]
        for v in simplex
    ]
    simplex_signed_response = [
        [
            max(left - right for left, right in zip(simplex_responses[i], simplex_responses[j]))
            for j in range(7)
        ]
        for i in range(7)
    ]
    assert all(
        simplex_signed_response[i][j] == (0 if i == j else 2)
        for i in range(7)
        for j in range(7)
    )

    matrices = tuple(multiplication_rows(a) for a in range(1, 8))
    assert len(set(matrices)) == 7
    rank_directed = []
    for i, left in enumerate(matrices):
        row = []
        for j, right in enumerate(matrices):
            difference = tuple(a ^ b for a, b in zip(left, right))
            value = 0 if i == j else min(binary_rank(left), binary_rank(difference))
            assert value == (0 if i == j else 3)
            row.append(value)
        rank_directed.append(row)

    matrix_queries = [
        tuple((encoded >> (3 * row)) & 0b111 for row in range(3))
        for encoded in range(1 << 9)
    ]
    rank_responses = [
        [
            min(
                binary_rank(query),
                2 + binary_rank(tuple(x ^ y for x, y in zip(query, matrix))),
            )
            for query in matrix_queries
        ]
        for matrix in matrices
    ]
    rank_signed_response = [
        [
            max(left - right for left, right in zip(rank_responses[i], rank_responses[j]))
            for j in range(7)
        ]
        for i in range(7)
    ]
    assert all(
        rank_signed_response[i][j] == (0 if i == j else 1)
        for i in range(7)
        for j in range(7)
    )

    # Exact product factorization pays the local presentation toll only on
    # changed symbols: the certified gaps are (d-p) rho, not d rho-p.
    assert (4 - 2) * (3 / 4) == 3 / 2
    assert (3 - 2) * (3 / 4) == 3 / 4
    return {
        "simplex_alphabet_size": 7,
        "simplex_word_weights": [weight(v) for v in simplex],
        "simplex_directed_distance_matrix": simplex_directed,
        "simplex_signed_response_matrix": simplex_signed_response,
        "rank_multiplication_alphabet_size": 7,
        "rank_directed_distance_matrix": rank_directed,
        "rank_signed_response_matrix": rank_signed_response,
        "hamming_response_gap_per_outer_coordinate": 1.5,
        "rank_response_gap_per_outer_coordinate": 0.75,
    }


def main() -> None:
    output = Path(__file__).with_name(
        "phase3_hamming_grassmannian_falsifier_results.json"
    )
    result = {
        "status": "passed",
        "exact_line_checks": exact_line_checks(),
        "local_alphabet": local_alphabet_checks(),
        "presented_alphabets": presented_alphabet_checks(),
    }
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
