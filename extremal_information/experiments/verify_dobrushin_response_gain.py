#!/usr/bin/env python3
"""Exact finite checks for the suffix-row response-gain identity."""

from fractions import Fraction
from itertools import product
from random import Random


def matmul(a, b):
    return tuple(tuple(sum(a[i][k] * b[k][j] for k in range(len(b)))
                       for j in range(len(b[0]))) for i in range(len(a)))


def identity(r):
    return tuple(tuple(Fraction(i == j) for j in range(r)) for i in range(r))


def tv(p, q):
    return sum(abs(x - y) for x, y in zip(p, q)) / 2


def hilbert(v):
    return (max(v) - min(v)) / 2


def mv(a, v):
    return tuple(sum(row[k] * v[k] for k in range(len(v))) for row in a)


def random_markov(rng, r, denominator=5):
    rows = []
    for _ in range(r):
        cuts = sorted([0, denominator] +
                      [rng.randrange(denominator + 1) for _ in range(r - 1)])
        weights = [cuts[i + 1] - cuts[i] for i in range(r)]
        rng.shuffle(weights)
        rows.append(tuple(Fraction(x, denominator) for x in weights))
    return tuple(rows)


def formula(mats):
    r = len(mats[0])
    suffixes = [None] * len(mats)
    cur = identity(r)
    for s in range(len(mats) - 1, -1, -1):
        suffixes[s] = cur
        cur = matmul(cur, mats[s])
    return max(sum(tv(q[i], q[j]) for q in suffixes)
               for i in range(r) for j in range(r))


def brute(mats):
    r = len(mats[0])
    vertices = [tuple(Fraction(2 * bit) for bit in bits)
                for bits in product((0, 1), repeat=r)]
    best = Fraction(0)
    for residuals in product(vertices, repeat=len(mats)):
        e = tuple(Fraction(0) for _ in range(r))
        for p, eta in zip(mats, residuals):
            e = tuple(x + y for x, y in zip(mv(p, e), eta))
        best = max(best, hilbert(e))
    return best


def main():
    rng = Random(20260816)
    checked = 0
    for r, depth, trials in ((2, 5, 120), (3, 3, 90)):
        for _ in range(trials):
            mats = tuple(random_markov(rng, r) for _ in range(depth))
            assert brute(mats) == formula(mats)
            checked += 1

    p = ((Fraction(3, 4), Fraction(1, 4)),
         (Fraction(1, 4), Fraction(3, 4)))
    for depth in range(1, 9):
        assert formula((p,) * depth) == sum(Fraction(1, 2) ** k
                                             for k in range(depth))
        checked += 1
    print(f"exact suffix-row gain checks: {checked}")


if __name__ == "__main__":
    main()
