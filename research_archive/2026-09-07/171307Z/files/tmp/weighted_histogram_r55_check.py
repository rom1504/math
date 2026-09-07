#!/usr/bin/env python3
"""Exact/numerical checks for the Wave 55 weighted-histogram memo."""

from fractions import Fraction
from itertools import combinations
from math import comb, sqrt

import numpy as np


def shell_count(n, m, j):
    k = n - m
    if not (0 <= j <= m and 0 <= m - j <= k):
        return 0
    return comb(m, j) * comb(k, m - j)


def d_core(n, m, ell):
    return comb(m, ell) * comb(n - ell, m - ell)


def lam(n, m, ell, degree):
    if degree > ell:
        return Fraction(0)
    ans = Fraction(1)
    for a in range(degree):
        ans *= Fraction((ell - a) * (n - m - a), (m - a) * (n - ell - a))
    return ans


def all_msets(n, m):
    return [frozenset(x) for x in combinations(range(n), m)]


def kernel_matrix(n, m, ell):
    sets = all_msets(n, m)
    den = d_core(n, m, ell)
    return np.array(
        [[comb(len(s & t), ell) / den for t in sets] for s in sets],
        dtype=float,
    )


def check_row_sum():
    print("row-sum identities")
    for n, m in [(12, 7), (14, 8), (16, 10)]:
        for ell in range(0, m):
            lhs = sum(
                shell_count(n, m, j) * comb(j, ell)
                for j in range(max(0, 2 * m - n), m + 1)
            )
            rhs = d_core(n, m, ell)
            assert lhs == rhs
        print(f"  n={n:2d}, m={m:2d}: all ell exact")


def check_spectrum():
    n, m = 8, 5
    weights = [(2, 0.4), (4, 0.6)]
    mat = sum(w * kernel_matrix(n, m, ell) for ell, w in weights)
    actual = np.linalg.eigvalsh(mat)[::-1]
    expected = []
    max_degree = min(m, n - m)
    for degree in range(max_degree + 1):
        eig = sum(w * float(lam(n, m, ell, degree)) for ell, w in weights)
        mult = comb(n, degree) - (comb(n, degree - 1) if degree else 0)
        expected.extend([eig] * mult)
    expected = np.array(sorted(expected, reverse=True))
    err = float(np.max(np.abs(actual - expected)))
    assert err < 2e-14
    assert np.max(np.abs(mat.sum(axis=1) - 1)) < 2e-14
    print("mixture spectrum")
    print(f"  n={n}, m={m}, weights={weights}, max error={err:.3e}")
    for degree in range(max_degree + 1):
        eig = sum(w * float(lam(n, m, ell, degree)) for ell, w in weights)
        print(f"  degree {degree}: Lambda={eig:.12f}")


def check_rearrangement():
    n, m = 7, 4
    weights = [(1, Fraction(1, 3)), (3, Fraction(2, 3))]
    sets = all_msets(n, m)
    base = sets[0]

    def q(j):
        return sum(w * Fraction(comb(j, ell), d_core(n, m, ell)) for ell, w in weights)

    vals = [q(len(base & t)) for t in sets if t != base]
    vals.sort(reverse=True)
    for r in [1, 2, 3]:
        brute = max(sum(x) for x in combinations(vals, r))
        predicted = sum(vals[:r])
        assert brute == predicted
    print("rearrangement capacity")
    print("  brute force agrees with filling intersection shells for r=1,2,3")


def family_retention(n, m, family, ell):
    den = d_core(n, m, ell)
    total = sum(comb(len(s & t), ell) for s in family for t in family)
    return Fraction(total, len(family) * den)


def check_convex_extraction():
    n, m = 10, 6
    family = [s for s in all_msets(n, m) if 0 in s]
    weights = [(2, Fraction(2, 5)), (4, Fraction(3, 5))]
    rows = []
    for ell, w in weights:
        p = family_retention(n, m, family, ell)
        l1, l2 = lam(n, m, ell, 1), lam(n, m, ell, 2)
        excess = p - l2
        den = (1 - l1) + n * (l1 - l2)
        rows.append((w, excess, den, excess / den))
    mixed_excess = sum(w * excess for w, excess, _, _ in rows)
    mixed_den = sum(w * den for w, _, den, _ in rows)
    mixed_ratio = mixed_excess / mixed_den
    assert min(x[3] for x in rows) <= mixed_ratio <= max(x[3] for x in rows)
    print("extraction convexity")
    for (ell, _), (_, excess, den, ratio) in zip(weights, rows):
        print(
            f"  ell={ell}: excess={float(excess):.9f}, "
            f"D={float(den):.9f}, ratio={float(ratio):.9f}"
        )
    print(f"  mixture ratio={float(mixed_ratio):.9f}")


def tuned_threshold_table():
    # p=3/5, alpha=9/25 gives beta_*=9/20 exactly.
    print("tuned one-threshold capacity")
    print("  n   ell  best_s  max B_s K_ell(s)  lambda2-h   ratio   sqrt(n)*max")
    for n in [100, 200, 400, 800, 1600, 3200]:
        assert n % 100 == 0
        m, ell = 3 * n // 5, 9 * n // 25
        den = d_core(n, m, ell)
        ball = 0
        best_val, best_s = Fraction(0), None
        for s in range(m - 1, ell - 1, -1):
            ball += shell_count(n, m, s)
            value = Fraction(ball * comb(s, ell), den)
            if value > best_val:
                best_val, best_s = value, s
        baseline = lam(n, m, ell, 2) - Fraction(comb(m, ell), den)
        assert baseline > 0
        print(
            f"  {n:4d} {ell:5d} {best_s:7d} {float(best_val):18.12f} "
            f"{float(baseline):11.9f} {float(best_val / baseline):7.4f} "
            f"{sqrt(n) * float(best_val):11.7f}"
        )


if __name__ == "__main__":
    check_row_sum()
    check_spectrum()
    check_rearrangement()
    check_convex_extraction()
    tuned_threshold_table()
