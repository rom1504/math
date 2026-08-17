#!/usr/bin/env python3
"""Finite exact checks for the PC.3 completed-parent limit."""

from fractions import Fraction
from itertools import product
from math import sqrt


def signs(word: str):
    return tuple(1 if c == "+" else -1 for c in word)


A = signs("+-----+--+-----+")
B = signs("+--+--------+--+")
C = signs("+++-+-++++-+-+++")


def mean_product(u, v):
    return Fraction(sum(x * y for x, y in zip(u, v)), len(u))


assert mean_product(A, B) == Fraction(1, 2)
assert mean_product(A, C) == 0
assert mean_product(B, C) == Fraction(-1, 2)

counts = {}
for a, b, c in zip(A, B, C):
    xy = (a * b, a * c)
    counts[xy] = counts.get(xy, 0) + 1
assert counts == {(1, 1): 4, (1, -1): 8, (-1, 1): 4}

local = tuple((xy, Fraction(k, 16)) for xy, k in counts.items())


def exact_L(j: int) -> Fraction:
    best = Fraction(0)
    # Enumerate endpoint coefficients and the product row law.  This is only
    # intended for small j; the proof gives the asymptotic statement.
    for eps0 in (-1, 1):
        for coeffs in product(product((-1, 1), repeat=2), repeat=j):
            value = Fraction(0)
            for rows in product(local, repeat=j):
                prob = Fraction(1)
                for (x, y), weight in rows:
                    prob *= weight
                total = eps0 + sum(
                    coeffs[t][0] * rows[t][0][0]
                    + coeffs[t][1] * rows[t][0][1]
                    for t in range(j)
                )
                value += prob * abs(total)
            best = max(best, value)
    return best


for j in range(1, 5):
    value = exact_L(j)
    lower = Fraction(j, 2) + 1
    upper = float(lower) + sqrt(11 * j / 4)
    assert value >= lower
    assert float(value) <= upper + 1e-12
    print(f"j={j}: L_j={value}, L_j/(2j+1)={float(value/(2*j+1)):.8f}")

print("PC.3 completed-parent limit checks: PASS")
