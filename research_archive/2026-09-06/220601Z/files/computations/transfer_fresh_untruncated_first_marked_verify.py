#!/usr/bin/env python3
"""Integer checks of untruncated first-marked identities and ascent constants.

No asymptotic or Gaussian conclusion is inferred from these finite checks.
Python standard library only; output is a short replay summary.
"""

from fractions import Fraction
from itertools import combinations, product
import random


def matvec(a, x):
    return [sum(t * s for t, s in zip(row, x)) for row in a]


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def verify(a, enumerate_cube):
    n = len(a)
    m = n - 1
    triples = list(combinations(range(n), 3))
    coefficients = []
    for i in range(n):
        ci = []
        for p, q, r in triples:
            ci.append(2 * (a[i][p] * a[p][q] * a[p][r]
                           + a[i][q] * a[p][q] * a[q][r]
                           + a[i][r] * a[p][r] * a[q][r]))
        coefficients.append(ci)
        assert sum(c * c for c in ci) <= 6 * m ** 3
        for seed in range(n):
            influence_numerator = sum(c * c for c, t in zip(ci, triples)
                                      if seed in t)
            assert influence_numerator <= 18 * m * (n - 2)

    if not enumerate_cube:
        return 0, 0

    spins = list(product((-1, 1), repeat=n))
    energies = []
    bilinear_cap = 0
    d_cov_sum = [[0] * n for _ in range(n)]
    for s in spins:
        field = matvec(a, s)
        energies.append(dot(s, field) // 2)
        bilinear_cap = max(bilinear_cap, sum(abs(t) for t in field))
        d_numerator = [s[i] * (field[i] ** 2 - m) for i in range(n)]
        y_numerator = matvec(a, d_numerator)
        for i in range(n):
            polynomial_value = sum(c * s[p] * s[q] * s[r]
                                   for c, (p, q, r)
                                   in zip(coefficients[i], triples))
            assert y_numerator[i] == polynomial_value
            for j in range(n):
                d_cov_sum[i][j] += d_numerator[i] * d_numerator[j]

    q_cap = max(abs(e) for e in energies)
    assert bilinear_cap <= 4 * q_cap
    for i in range(n):
        for j in range(n):
            a_squared = sum(a[i][k] * a[k][j] for k in range(n))
            expected = 4 * a_squared + (2 * m * (m - 3) if i == j else 0)
            assert d_cov_sum[i][j] == len(spins) * expected

    ascent_count = 0
    for s in spins:
        field = matvec(a, s)
        for orientation in (-1, 1):
            unstable = [orientation * s[i] * field[i] < 0 for i in range(n)]
            w = sum(-orientation * s[i] * field[i]
                    for i in range(n) if unstable[i])
            v = [s[i] if unstable[i] else 0 for i in range(n)]
            p = Fraction(w, 4 * q_cap)
            assert 0 <= p <= 1
            actual_gain = 2 * p * w + 2 * orientation * p * p * dot(v, matvec(a, v))
            assert actual_gain >= Fraction(w * w, 4 * q_cap)
            ascent_count += 1
    return len(spins), ascent_count


def main():
    rng = random.Random(20260906)
    cases = cube_states = ascent_count = 0
    for n in range(2, 17):
        for trial in range(8):
            a = [[0] * n for _ in range(n)]
            for i in range(n):
                for j in range(i):
                    a[i][j] = a[j][i] = 1 if trial == 0 else rng.choice((-1, 1))
            states, ascents = verify(a, enumerate_cube=(n <= 8))
            cube_states += states
            ascent_count += ascents
            cases += 1
    print("PASS: {} signings, orders 2..16; {} exact cube states; {} oriented ascent checks."
          .format(cases, cube_states, ascent_count))
    print("Verified cubic coefficients, covariance D, variance/influence bounds, and cap-only ascent.")


if __name__ == "__main__":
    main()
