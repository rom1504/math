#!/usr/bin/env python3
"""Exact restriction/jackknife audit for the certified order-nine minimizer."""

from fractions import Fraction
from itertools import combinations, product


A = (
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

# Doubled exact minima, from the certified small-order data in the project.
QMIN = {1: 0, 2: 2, 3: 6, 4: 8, 5: 8, 6: 10, 7: 18, 8: 20, 9: 24}


def q_of_subset(vertices):
    vertices = tuple(vertices)
    if len(vertices) <= 1:
        return 0
    # Fix the first spin to +1 (global sign invariance).
    best = 0
    for tail in product((-1, 1), repeat=len(vertices) - 1):
        x = (1,) + tail
        energy = 0
        for p in range(len(vertices)):
            for q in range(p + 1, len(vertices)):
                energy += 2 * A[vertices[p]][vertices[q]] * x[p] * x[q]
        best = max(best, abs(energy))
    return best


def beta(k, m):
    return Fraction(m * (m - 1), k * (k - 1))


def main():
    by_order = {}
    qmap = {}
    expected = {
        1: (0, 0, Fraction(0)),
        2: (2, 2, Fraction(2)),
        3: (6, 6, Fraction(6)),
        4: (8, 12, Fraction(184, 21)),
        5: (8, 16, Fraction(248, 21)),
        6: (10, 22, Fraction(118, 7)),
        7: (18, 22, Fraction(21)),
        8: (24, 24, Fraction(24)),
        9: (24, 24, Fraction(24)),
    }
    universe = range(9)
    for m in range(1, 10):
        subsets = list(combinations(universe, m))
        values = [q_of_subset(s) for s in subsets]
        qmap.update(zip(subsets, values))
        by_order[m] = values
        avg = Fraction(sum(values), len(values))
        assert (min(values), max(values), avg) == expected[m]
        excess = avg - QMIN[m]
        residual = avg - beta(9, m) * QMIN[9] if m >= 2 else avg
        print(
            f"m={m}: count={len(values)}, min={min(values)}, max={max(values)}, "
            f"F={avg}, eps={excess}, J_9m={residual}"
        )

    # Verify the exact two-level cocycle for every 2 <= m < k <= 9.
    for k in range(3, 10):
        Fk = Fraction(sum(by_order[k]), len(by_order[k]))
        J9k = Fk - beta(9, k) * QMIN[9]
        for m in range(2, k):
            Fm = Fraction(sum(by_order[m]), len(by_order[m]))
            J9m = Fm - beta(9, m) * QMIN[9]
            averaged_local = Fm - beta(k, m) * Fk
            assert J9m == averaged_local + beta(k, m) * J9k
    print("all two-level cocycles: PASS")

    # Conditional 8 -> 7 fresh-ground averages and exact variances.
    print("conditional 8->7 data:")
    for U in combinations(universe, 8):
        vals = [qmap[tuple(v for v in U if v != i)] for i in U]
        mean = Fraction(sum(vals), len(vals))
        var = sum((Fraction(v) - mean) ** 2 for v in vals) / len(vals)
        print(f"  missing={next(i for i in universe if i not in U)} mean={mean} var={var}")

    flat_u = tuple(i for i in universe if i != 7)
    assert qmap[flat_u] == 24
    assert {qmap[tuple(v for v in flat_u if v != i)] for i in flat_u} == {22}
    print("flat order-8 child restriction profile:")
    for m in range(2, 9):
        vals = [qmap[S] for S in combinations(flat_u, m)]
        mean = Fraction(sum(vals), len(vals))
        var = sum((Fraction(v) - mean) ** 2 for v in vals) / len(vals)
        print(
            f"  m={m} min={min(vals)} max={max(vals)} mean={mean} "
            f"eps={mean-QMIN[m]} var={var}"
        )


if __name__ == "__main__":
    main()
