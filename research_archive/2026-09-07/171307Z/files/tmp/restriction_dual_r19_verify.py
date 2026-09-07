#!/usr/bin/env python3
"""Exact checks for the Wave 19 restriction-game memo; standard library only."""

from fractions import Fraction
from itertools import combinations, product
import random


A9 = (
    (0, 1, -1, 1, 1, 1, 1, -1, -1),
    (1, 0, 1, -1, -1, -1, 1, -1, -1),
    (-1, 1, 0, 1, 1, 1, 1, 1, -1),
    (1, -1, 1, 0, 1, 1, 1, -1, 1),
    (1, -1, 1, 1, 0, 1, -1, 1, -1),
    (1, -1, 1, 1, 1, 0, -1, -1, -1),
    (1, 1, 1, 1, -1, -1, 0, 1, 1),
    (-1, -1, 1, -1, 1, -1, 1, 0, -1),
    (-1, -1, -1, 1, -1, -1, 1, -1, 0),
)


def energy_on(a, x, vertices):
    return 2 * sum(
        a[i][j] * x[i] * x[j] for i, j in combinations(vertices, 2)
    )


def a9_lp_certificate():
    n = 9
    rows = [tuple(j for j in range(n) if j != i) for i in range(n)]
    columns = []
    metadata = []
    for tail in product((-1, 1), repeat=n - 1):
        x = (1,) + tail
        for sigma in (1, -1):
            columns.append(tuple(sigma * energy_on(A9, x, s) for s in rows))
            metadata.append((x, sigma))

    assert [max(col[i] for col in columns) for i in range(n)] == [24] * n

    pi_num = (4, 2, 4, 0, 4, 4, 2, 5, 0)
    primal_values = [
        Fraction(sum(pi_num[i] * col[i] for i in range(n)), 25)
        for col in columns
    ]
    assert max(primal_values) == Fraction(464, 25)

    dual_atoms = (
        ((1, -1, 1, -1, -1, -1, 1, -1, -1), -1, 388),
        ((1, -1, 1, -1, 1, -1, 1, -1, 1), -1, 64),
        ((1, -1, 1, 1, 1, 1, -1, -1, -1), 1, 149),
        ((1, -1, 1, 1, 1, 1, -1, 1, -1), 1, 173),
        ((1, 1, 1, -1, -1, 1, -1, 1, 1), -1, 160),
        ((1, 1, 1, 1, -1, -1, -1, 1, 1), -1, 340),
        ((1, 1, 1, 1, 1, 1, 1, -1, -1), 1, 51),
        ((1, 1, 1, 1, 1, 1, 1, 1, -1), 1, 75),
    )
    assert sum(weight for _, _, weight in dual_atoms) == 1400
    dual_rows = []
    for vertices in rows:
        numerator = sum(
            weight * sigma * energy_on(A9, x, vertices)
            for x, sigma, weight in dual_atoms
        )
        dual_rows.append(Fraction(numerator, 1400))
    assert dual_rows == [
        Fraction(464, 25), Fraction(464, 25), Fraction(464, 25),
        Fraction(488, 25), Fraction(464, 25), Fraction(464, 25),
        Fraction(464, 25), Fraction(464, 25), Fraction(464, 25),
    ]
    print("A9 LP: V_ad=24, V_hid=464/25; rational primal/dual certificates pass")


def falling(a, k):
    out = 1
    for j in range(k):
        out *= a - j
    return out


def variance_audit():
    rng = random.Random(1901)
    for n in range(4, 9):
        for m in range(2, n):
            subsets = list(combinations(range(n), m))
            p = {k: Fraction(falling(m, k), falling(n, k)) for k in range(2, 5)}
            for _ in range(12):
                a = [[0] * n for _ in range(n)]
                for i, j in combinations(range(n), 2):
                    a[i][j] = a[j][i] = rng.choice((-1, 1))
                x = [rng.choice((-1, 1)) for _ in range(n)]
                w = {(i, j): a[i][j] * x[i] * x[j] for i, j in combinations(range(n), 2)}
                ys = [sum(w[e] for e in combinations(s, 2)) for s in subsets]
                mean = sum(map(Fraction, ys)) / len(ys)
                var = sum((Fraction(y) - mean) ** 2 for y in ys) / len(ys)
                W = sum(w.values())
                r = [sum(a[i][j] * x[i] * x[j] for j in range(n) if j != i) for i in range(n)]
                N = n * (n - 1) // 2
                rhs = (
                    (p[3] - p[4]) * sum(z * z for z in r)
                    + N * (p[2] - 2 * p[3] + p[4])
                    + (p[4] - p[2] * p[2]) * W * W
                )
                assert mean == p[2] * W
                assert var == rhs
    print("selector variance: exact exhaustive subset audits pass")


if __name__ == "__main__":
    a9_lp_certificate()
    variance_audit()
