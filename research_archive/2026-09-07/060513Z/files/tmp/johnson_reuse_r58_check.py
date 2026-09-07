#!/usr/bin/env python3
"""Exact audit of the Wave 58 common-state Johnson-reuse wall."""

from fractions import Fraction
from itertools import combinations, product
from math import comb

import numpy as np


C = np.asarray(
    [
        [0, 1, 1, 1, 1, 1],
        [1, 0, -1, -1, 1, 1],
        [1, -1, 0, 1, -1, 1],
        [1, -1, 1, 0, 1, -1],
        [1, 1, -1, 1, 0, -1],
        [1, 1, 1, -1, -1, 0],
    ],
    dtype=np.int64,
)
Z0 = np.asarray((1, -1, 1, 1, -1, 1), dtype=np.int64)
Z1 = np.asarray((1, -1, -1, 1, 1, -1), dtype=np.int64)
E = ((0, 1), (0, 2), (0, 3), (0, 4), (0, 5))


def projective_spins(n: int):
    for tail in product((-1, 1), repeat=n - 1):
        yield np.asarray((1,) + tail, dtype=np.int64)


def child_sector_audit() -> None:
    values = [int(x @ C @ x) for x in projective_spins(6)]
    q = max(map(abs, values))
    one = np.ones(6, dtype=np.int64)
    assert q == 10 and int(one @ C @ one) == q
    assert tuple(C @ one) == (5, 1, 1, 1, 1, 1)

    # The per-vertex selections all give the displayed star.
    fields = C @ one
    selected = set()
    for i in range(6):
        positive = [(min(i, j), max(i, j)) for j in range(6)
                    if i != j and C[i, j] == 1]
        if i == 0:
            choice = positive
        else:
            choice = [(0, i)]
        assert len(choice) == int(fields[i])
        assert all(C[a, b] == 1 for a, b in choice)
        selected.update(choice)
    assert tuple(sorted(selected)) == E
    assert len(E) == q // 2

    records = []
    for z in (Z0, Z1):
        disagreement = [(i, j) for i, j in combinations(range(6), 2)
                        if z[i] * z[j] == -1]
        cut_sum = sum(int(C[i, j]) for i, j in disagreement)
        block_sum = sum(int(C[i, j] * z[i] * z[j]) for i, j in E)
        local_parent_energy = int(z @ C @ z)
        delta = q - local_parent_energy
        assert delta == 4 * cut_sum
        assert block_sum == len(E) - 2 * sum(
            z[i] * z[j] == -1 for i, j in E
        )
        assert 2 * block_sum <= len(E)
        records.append((cut_sum, block_sum, delta))
    assert records == [(0, 1, 0), (1, -1, 4)]
    print("child templates: cap=10, legal |E|=5, records", records)


def johnson_audit() -> None:
    n, m, ell = 11, 6, 3
    omega = [frozenset(s) for s in combinations(range(n), m)]
    family = [
        frozenset(s) for s in (
            (0, 1, 3, 4, 6, 7),
            (0, 1, 2, 6, 7, 8),
            (0, 3, 4, 5, 7, 10),
            (0, 4, 5, 6, 8, 10),
            (2, 3, 5, 6, 7, 9),
            (1, 5, 7, 8, 9, 10),
            (0, 1, 2, 3, 4, 5),
            (1, 3, 6, 8, 9, 10),
        )
    ]
    assert all(s in omega for s in family)
    N, r = len(omega), len(family)
    d = comb(m, ell) * comb(n - ell, m - ell)
    h = Fraction(1, comb(n - ell, m - ell))
    lam2 = Fraction(
        comb(m - 2, ell - 2) * comb(n - ell - 2, m - ell), d
    )
    weighted = sum(comb(len(s & t), ell) for s in family for t in family)
    P = Fraction(weighted, r * d)
    assert (N, d, h, lam2, P) == (
        462, 1120, Fraction(1, 56), Fraction(1, 14), Fraction(3, 112)
    )
    assert P < lam2

    rows = []
    for threshold in range(ell, m):
        D = sum(comb(m, j) * comb(n - m, m - j)
                for j in range(threshold, m))
        R = Fraction(D * comb(threshold, ell), d)
        pairs = sum(
            s != t and len(s & t) >= threshold
            for s in family for t in family
        )
        Pi = Fraction(pairs, r * D)
        assert P >= h + R * Pi
        rows.append((R, threshold, D, pairs, Pi))
    R, s_star, D, pairs, Pi = max(rows)
    critical = (lam2 - h) / R
    assert (s_star, D, R, pairs, Pi, critical) == (
        4, 180, Fraction(9, 14), 12, Fraction(1, 120), Fraction(1, 12)
    )
    assert Pi < critical

    # Exact random-r-subset expectations from ordered-pair counting.
    expected_P = h + (1 - h) * Fraction(r - 1, N - 1)
    expected_P_direct = Fraction(
        r * comb(m, ell)
        + Fraction(r * (r - 1), N * (N - 1))
        * sum(
            comb(len(s & t), ell)
            for s in omega for t in omega if s != t
        ),
        r * d,
    )
    assert expected_P == expected_P_direct
    expected_Pi = Fraction(r - 1, N - 1)
    assert expected_Pi == Fraction(7, 461)
    print(
        "Johnson audit: P3=3/112 < lambda2=1/14; "
        "s*=4, Pi=1/120 < critical=1/12"
    )


if __name__ == "__main__":
    child_sector_audit()
    johnson_audit()
    print("PASS exact child-sector compatibility and full Johnson histogram")
