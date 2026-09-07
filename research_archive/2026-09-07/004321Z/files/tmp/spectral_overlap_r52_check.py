#!/usr/bin/env python3
"""Exact checks for the Wave 52 spectral-overlap audit.

This script verifies the core-load formula, the exact level-two surplus,
the direct/spectral hybrid degree bounds, and a small exact-minimizer atlas.
"""

from __future__ import annotations

import itertools
import math
import sys
from collections import Counter
from fractions import Fraction

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from box_triple_r50_check import A as A10
from coarea_alignment_r47_check import exact_boundary_data
from envelope_block_cover_r27 import A6, A8, A9


def falling(x: int, j: int) -> int:
    ans = 1
    for i in range(j):
        ans *= x - i
    return ans


def eigenvalue(n: int, m: int, ell: int, j: int) -> Fraction:
    if ell < j or n - m < j:
        return Fraction(0)
    return Fraction(
        falling(ell, j) * falling(n - m, j),
        falling(m, j) * falling(n - ell, j),
    )


def family_core_data(n: int, m: int, family, ell: int):
    states = [frozenset(s) for s in family]
    r = len(states)
    loads = Counter(
        core
        for s in states
        for core in itertools.combinations(sorted(s), ell)
    )
    d = math.comb(m, ell) * math.comb(n - ell, m - ell)
    p = Fraction(sum(c * c for c in loads.values()), r * d)
    rho = Fraction(math.comb(m, ell), math.comb(n, ell))
    a = Fraction(r, math.comb(n, m))
    # p = a E_R (mu_R/rho)^2, including zero-load cores.
    normalized_second = Fraction(
        sum(c * c for c in loads.values()) * math.comb(n, ell),
        r * r * math.comb(m, ell) ** 2,
    )
    assert p == a * normalized_second
    assert a <= p <= a / rho
    return p, a, rho, loads


def symbolic_identities() -> None:
    checked = 0
    for n in range(7, 17):
        for m in range(3, n - 1):
            for ell in range(2, m):
                d = math.comb(m, ell) * math.comb(n - ell, m - ell)
                lam2 = eigenvalue(n, m, ell, 2)
                baseline = (
                    math.comb(m - 2, ell - 2)
                    * math.comb(n - ell - 2, m - ell)
                )
                assert lam2 * d == baseline
                checked += 1
    print(f"symbolic lambda2*d_ell identity: {checked} parameter triples")


def exhaustive_family_identities() -> None:
    n, m = 5, 3
    universe = [frozenset(s) for s in itertools.combinations(range(n), m)]
    checked = 0
    for mask in range(1, 1 << len(universe)):
        family = [universe[i] for i in range(len(universe)) if mask >> i & 1]
        for ell in (1, 2):
            family_core_data(n, m, family, ell)
            checked += 1
    print(f"exhaustive core-load formula/bounds: {checked} family-scales")


def aggregate_retention(f, rows, selectors, cap: int, ell: int):
    n = max(max(s) for s in selectors) + 1
    m = len(selectors[0])
    d = math.comb(m, ell) * math.comb(n - ell, m - ell)
    numerator = 0
    denominator = 0
    active_sizes = []
    for iz in range(len(rows)):
        ids = np.flatnonzero(f[iz])
        r = len(ids)
        if rows[iz] > cap or not r:
            continue
        active_sizes.append(r)
        moment = sum(
            math.comb(len(set(selectors[i]) & set(selectors[j])), ell)
            for i in ids for j in ids
        )
        numerator += r * moment
        denominator += r * r
    assert denominator
    return Fraction(numerator, d * denominator), active_sizes


def finite_atlas() -> None:
    tested = 0
    failures = []
    high_density = []
    for name, a in (("A6", A6), ("A8", A8), ("A9", A9), ("A10", A10)):
        n = len(a)
        m_values = range(3, n) if name != "A10" else (6,)
        for m in m_values:
            z, selectors, q, _, rows, _, _, f, degree = exact_boundary_data(a, m)
            active = [i for i, az in enumerate(degree) if az > 0]
            cap = min(int(rows[i]) for i in active)
            p2, sizes = aggregate_retention(f, rows, selectors, cap, 2)
            lam2 = eigenvalue(n, m, 2, 2)
            passed = p2 >= lam2
            tested += 1
            if not passed:
                failures.append((name, n, m, cap, sizes, p2, lam2))
            if 2 * m > n and n - m >= 2:
                self_loop = Fraction(1, math.comb(n - 2, m - 2))
                high_density.append(
                    (name, n, m, cap, sizes, p2, lam2, self_loop, passed)
                )
    assert failures == [
        ("A8", 8, 3, 8, [8, 8], Fraction(1, 6), Fraction(2, 9)),
        ("A9", 9, 3, 16, [14, 14], Fraction(11, 49), Fraction(5, 21)),
    ]
    assert all(row[-1] and row[-2] >= row[-3] for row in high_density)
    print(
        f"finite sanity atlas (self-loop contaminated, not signing evidence): "
        f"{tested} cases; failures={failures}"
    )
    print("m>n/2 rows (all have h2>=lambda2 before ground structure):")
    for row in high_density:
        print("  " + " | ".join(map(str, row)))


def order_ten_pair_surplus() -> None:
    n, m = 10, 6
    z, selectors, q, _, rows, _, _, f, degree = exact_boundary_data(A10, m)
    active = [i for i, az in enumerate(degree) if az > 0]
    cap = min(int(rows[i]) for i in active)
    low = [i for i in active if int(rows[i]) == cap]
    assert cap == 10 and len(low) == 2
    families = [tuple(np.flatnonzero(f[i])) for i in low]
    assert families[0] == families[1] and len(families[0]) == 5
    family = [frozenset(selectors[i]) for i in families[0]]
    p2, a, rho2, loads = family_core_data(n, m, family, 2)
    lam2 = eigenvalue(n, m, 2, 2)
    d2 = math.comb(m, 2) * math.comb(n - 2, m - 2)
    baseline = math.comb(n - 4, m - 2)
    sum_squares = sum(c * c for c in loads.values())
    pair_surplus = sum_squares - baseline * len(family)
    assert p2 == Fraction(29, 1050)
    assert lam2 == Fraction(1, 70)
    assert pair_surplus == 70
    assert p2 - lam2 == Fraction(pair_surplus, len(family) * d2)
    intersection_histogram = Counter(
        len(s & t) for s in family for t in family if s != t
    )
    # In particular no favorable selectors are adjacent in the Johnson graph.
    assert intersection_histogram == Counter({4: 10, 2: 10})
    assert all(j != m - 1 for j in intersection_histogram)
    load_histogram = Counter(loads.values())
    print(
        "A10,m=6,row=10: "
        f"a={a}, rho2={rho2}, P2={p2}, lambda2={lam2}, "
        f"pair_surplus={pair_surplus}, load_hist={dict(sorted(load_histogram.items()))}, "
        f"distinct-overlap-hist={dict(sorted(intersection_histogram.items()))}"
    )


def hybrid_bound_examples() -> None:
    # If P >= lambda2, the direct core-load bound already gives
    # M >= rho_2 lambda_2, even when the spectral excess is exactly zero.
    for n in (40, 80, 160, 320):
        m = 3 * n // 5
        rho2 = Fraction(falling(m, 2), falling(n, 2))
        lam2 = eigenvalue(n, m, 2, 2)
        direct = rho2 * lam2
        claimed = Fraction(
            2 * (n - m) * (n - m - 1),
            n * (n - 1) * (n - 2) * (n - 3),
        )
        assert direct == claimed
        print(f"level-two direct threshold n={n},m={m}: rho2*lambda2={direct}")


def main() -> None:
    symbolic_identities()
    exhaustive_family_identities()
    order_ten_pair_surplus()
    finite_atlas()
    hybrid_bound_examples()
    print("PASS: exact core-load, surplus, regime, and finite-atlas checks")


if __name__ == "__main__":
    main()
