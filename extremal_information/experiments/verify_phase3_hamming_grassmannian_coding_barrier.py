#!/usr/bin/env python3
"""Exhaustive finite checks for the Hamming Grassmannian coding barrier."""

from itertools import combinations, product
from math import comb


def wt(x: int) -> int:
    return bin(x).count("1")


def line_distance(v: int, w: int) -> int:
    c = wt(v ^ w)
    return max(min(wt(v), c), min(wt(w), c))


def line_formula(v: int, w: int) -> int:
    lv = (0, v)
    lw = (0, w)
    directed_vw = max(min(wt(x ^ y) for y in lw) for x in lv)
    directed_wv = max(min(wt(y ^ x) for x in lv) for y in lw)
    return max(directed_vw, directed_wv)


def maximum_clique_size(vertices, adjacent) -> int:
    """Small exact Bron--Kerbosch maximum-clique search with bitsets."""
    n = len(vertices)
    nbr = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if adjacent(vertices[i], vertices[j]):
                nbr[i] |= 1 << j
                nbr[j] |= 1 << i

    best = 0

    def expand(size: int, candidates: int) -> None:
        nonlocal best
        if size + wt(candidates) <= best:
            return
        while candidates:
            if size + wt(candidates) <= best:
                return
            bit = candidates & -candidates
            i = bit.bit_length() - 1
            candidates ^= bit
            expand(size + 1, candidates & nbr[i])
        best = max(best, size)

    expand(0, (1 << n) - 1)
    return best


def binary_code_number(D: int, d: int) -> int:
    words = list(range(1 << D))
    return maximum_clique_size(words, lambda x, y: wt(x ^ y) >= d)


def line_packing_number(D: int, t: int) -> int:
    words = list(range(1, 1 << D))
    return maximum_clique_size(words, lambda x, y: line_distance(x, y) > t)


def intersection_volume(D: int, t: int, b: int) -> int:
    ans = 0
    for i in range(b + 1):
        for j in range(D - b + 1):
            if i + j <= t and b - i + j <= t:
                ans += comb(b, i) * comb(D - b, j)
    return ans


def predicted_line_ball(D: int, t: int, w: int) -> int:
    volume = sum(comb(D, i) for i in range(t + 1))
    b = wt(w)
    if b > t:
        return volume
    return 2 * volume - intersection_volume(D, t, b) - 1


def span_from_rows(rows):
    points = {0}
    for row in rows:
        points |= {x ^ row for x in tuple(points)}
    return frozenset(points)


def graph_subspace(k: int, columns):
    rows = []
    for i in range(k):
        suffix = sum(((col >> i) & 1) << j for j, col in enumerate(columns))
        rows.append((1 << i) | (suffix << k))
    return span_from_rows(rows)


def hausdorff(C, E) -> int:
    directed_ce = max(min(wt(c ^ e) for e in E) for c in C)
    directed_ec = max(min(wt(e ^ c) for c in C) for e in E)
    return max(directed_ce, directed_ec)


def verify_lines() -> None:
    for D in range(1, 6):
        nonzero = range(1, 1 << D)
        for v, w in combinations(nonzero, 2):
            assert line_distance(v, w) == line_formula(v, w)
        for t in range(D):
            A = binary_code_number(D, t + 1)
            P = line_packing_number(D, t)
            assert A - 1 <= P <= A, (D, t, A, P)
            for w in nonzero:
                actual = sum(line_distance(v, w) <= t for v in nonzero)
                assert actual == predicted_line_ball(D, t, w)
    print("line formula, balls, and A_2-1 <= P <= A_2: exhaustive D <= 5")


def verify_systematic_charts() -> None:
    cases = 0
    for k, L in ((1, 4), (2, 3), (3, 2)):
        alphabet = range(1 << k)
        matrices = list(product(alphabet, repeat=L))
        carriers = [graph_subspace(k, cols) for cols in matrices]
        for i, j in combinations(range(len(matrices)), 2):
            column_distance = sum(a != b for a, b in zip(matrices[i], matrices[j]))
            assert hausdorff(carriers[i], carriers[j]) <= column_distance
            cases += 1
    print(f"systematic-chart inequality: {cases:,} exhaustive matrix pairs")


if __name__ == "__main__":
    verify_lines()
    verify_systematic_charts()
