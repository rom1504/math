#!/usr/bin/env python3
"""Independent finite checks for the exact identities in center_power_r32.md."""

from fractions import Fraction
from itertools import combinations, product
from math import comb, log
import random


def mat_vec(a, x):
    return [sum(a[i][j] * x[j] for j in range(len(x))) for i in range(len(x))]


def row_square(a, x):
    return sum(v * v for v in mat_vec(a, x))


def mat_square(a):
    n = len(a)
    return [[sum(a[i][k] * a[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def conditional_checks():
    rng = random.Random(3202)
    n, m = 7, 4
    a = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            a[i][j] = a[j][i] = rng.choice((-1, 1))
    q = mat_square(a)
    assert all(q[i][i] == n - 1 for i in range(n))
    y = (1, -1, 1, 1)

    vals = []
    for outside in product((-1, 1), repeat=n - m):
        vals.append(row_square(a, y + outside))
    brute = Fraction(sum(vals), len(vals))
    formula = n * (n - 1) + sum(
        q[i][j] * y[i] * y[j]
        for i in range(m)
        for j in range(m)
        if i != j
    )
    assert brute == formula

    for h in range(m + 1):
        sphere_vals = []
        for flip in combinations(range(m), h):
            yy = list(y)
            for i in flip:
                yy[i] *= -1
            for outside in product((-1, 1), repeat=n - m):
                sphere_vals.append(row_square(a, tuple(yy) + outside))
        brute_h = Fraction(sum(sphere_vals), len(sphere_vals))
        if m > 1:
            theta = Fraction(m * (m - 1) - 4 * h * (m - h), m * (m - 1))
        else:
            theta = Fraction(1)
        predicted = n * (n - 1) + theta * (formula - n * (n - 1))
        assert brute_h == predicted


def cylinder_distance_check():
    # Two projective labels on S; outside completion bits are unrestricted.
    n = 6
    s_idx = (0, 2, 4)
    labels = ((1, 1, -1), (1, -1, -1))
    z = (1, -1, -1, 1, 1, -1)

    def ham(u, v):
        return sum(a != b for a, b in zip(u, v))

    restricted = tuple(z[i] for i in s_idx)
    lhs = min(min(ham(restricted, y), ham(restricted, tuple(-v for v in y))) for y in labels)
    full_best = n
    outside_idx = tuple(i for i in range(n) if i not in s_idx)
    for y in labels:
        for sign in (-1, 1):
            oriented = tuple(sign * v for v in y)
            for outside in product((-1, 1), repeat=len(outside_idx)):
                x = [0] * n
                for i, v in zip(s_idx, oriented):
                    x[i] = v
                for i, v in zip(outside_idx, outside):
                    x[i] = v
                full_best = min(full_best, ham(z, x), ham(z, tuple(-v for v in x)))
    assert lhs == full_best


def hamming_and_exponent_checks():
    # Exact volume bound on a representative finite pair.
    m, d = 120, 7
    volume = sum(comb(m, j) for j in range(d + 1))
    assert log(volume) <= d * log(2.718281828459045 * m / d)

    # Project exponent inequalities for several admissible (c0,eta).
    for c0, eta in ((0.02, 0.0), (0.10, 0.05), (0.24, 0.239)):
        assert 0 <= eta < c0 < 0.25
        d_exp = 0.5 - 2 * c0 + eta
        l_exp = 0.75 - c0
        tl_exp = 0.75 - c0 + eta
        ratio_exp = 0.25 + c0 - eta
        assert d_exp < l_exp
        assert tl_exp < 1
        assert ratio_exp > 0.25
        assert abs((1 + (0.25 + c0 - eta)) - (1.25 + c0 - eta)) < 1e-12

    # max_{0<=p<=1} p^(3/2)-p^2 = 27/256 at p=9/16.
    p = Fraction(9, 16)
    # sqrt(p)=3/4, so evaluate without floating-point powers.
    value = p * Fraction(3, 4) - p * p
    assert value == Fraction(27, 256)


def clique_identity_check():
    for b in range(3, 10):
        for z in product((-1, 1), repeat=b):
            magnetization = sum(z)
            # P_B z has coordinate M-z_i on B.
            brute = sum((magnetization - z_i) ** 2 for z_i in z)
            formula = (b - 2) * magnetization**2 + b
            assert brute == formula


if __name__ == "__main__":
    conditional_checks()
    cylinder_distance_check()
    hamming_and_exponent_checks()
    clique_identity_check()
    print("PASS center_power_r32: cylinder, noise, volume, exponents, clique identity")
