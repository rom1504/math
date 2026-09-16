#!/usr/bin/env python3
"""Exact audit of (10.442) on the strict denominator-672 witness.

All cut quantities are undoubled.  This is scratch validation only.
"""

from fractions import Fraction as F
from itertools import combinations, product


S = tuple(range(4))
T = tuple(range(4, 9))
V = S + T


def weights():
    w = {}
    # S is a positive K4.
    for i, j in combinations(S, 2):
        w[i, j] = F(9, 112)
    # All four cross rows are identical.
    cross_row = (F(1, 24), F(1, 24), F(1, 24), F(55, 672), F(29, 672))
    for i in S:
        for j, q in zip(T, cross_row):
            w[i, j] = q
    # T: a positive triangle, six equal negative spokes, and one positive edge.
    for i, j in combinations((4, 5, 6), 2):
        w[i, j] = F(83, 672)
    for i in (4, 5, 6):
        for j in (7, 8):
            w[i, j] = F(-55, 672)
    w[7, 8] = F(113, 672)
    assert len(w) == 36
    return w


W = weights()


def cut(ids, spin):
    pos = dict(zip(ids, spin))
    return sum(q for (i, j), q in W.items()
               if i in pos and j in pos and pos[i] != pos[j])


def total(ids):
    ids = set(ids)
    return sum(q for (i, j), q in W.items() if i in ids and j in ids)


def states(ids):
    return [(1,) + tail for tail in product((-1, 1), repeat=len(ids) - 1)]


def cross_bilinear(x, y):
    return sum(W[i, j] * x[u] * y[v]
               for u, i in enumerate(S) for v, j in enumerate(T))


def main():
    full = states(V)
    fvals = [(x, cut(V, x)) for x in full]
    endpoint = (1,) * len(S) + (-1,) * len(T)
    assert cut(V, endpoint) == 1
    assert min(q for _, q in fvals) == 0
    assert max(q for _, q in fvals) == 1
    assert [x for x, q in fvals if q == 0] == [(1,) * 9]
    assert [x for x, q in fvals if q == 1] == [endpoint]
    nonend = [q for x, q in fvals if x not in ((1,) * 9, endpoint)]
    assert min(nonend) == F(1, 112)
    assert max(nonend) == F(111, 112)

    child = []
    maxima = []
    minima = []
    for ids in (S, T):
        vals = [(x, cut(ids, x)) for x in states(ids)]
        lo = min(q for _, q in vals)
        hi = max(q for _, q in vals)
        child.append((total(ids), lo, hi))
        minima.append([x for x, q in vals if q == lo])
        maxima.append([x for x, q in vals if q == hi])

    aS, mS, MS = child[0]
    aT, mT, MT = child[1]
    assert (aS, mS, MS) == (F(27, 56), F(0), F(9, 28))
    assert (aT, mT, MT) == (F(1, 21), F(-55, 112), F(19, 112))
    assert MS < aS - mS and MT < aT - mT

    expected_x = [
        (1, -1, -1, 1),
        (1, -1, 1, -1),
        (1, 1, -1, -1),
    ]
    expected_y = [
        (1, -1, -1, -1, 1),
        (1, -1, -1, 1, -1),
        (1, -1, 1, -1, 1),
        (1, -1, 1, 1, -1),
        (1, 1, -1, -1, 1),
        (1, 1, -1, 1, -1),
    ]
    assert maxima == [expected_x, expected_y]
    zvals = {(x, y): cross_bilinear(x, y)
             for x in maxima[0] for y in maxima[1]}
    assert set(zvals.values()) == {F(0)}

    target = aS + aT - mS - mT
    lhs = 2 * (MS + MT) + max(abs(z) for z in zvals.values())
    assert target == F(49, 48)
    assert lhs == F(55, 56)
    assert target - lhs == F(13, 336)

    # Audit the corresponding doubled-energy comparison in the ledger.
    favorable_parent_energy = 2 * (aS + aT) - 4 * (MS + MT)
    child_min_energy = 2 * (mS + mT)
    assert favorable_parent_energy == F(-19, 21)
    assert child_min_energy == F(-55, 56)
    assert favorable_parent_energy - child_min_energy == F(13, 168)

    print('parent c/min_nonendpoint/max_nonendpoint', 1, min(nonend), max(nonend))
    print('child S (a,m,M)', child[0], 'minimizers', minima[0], 'maximizers', maxima[0])
    print('child T (a,m,M)', child[1], 'minimizers', minima[1], 'maximizers', maxima[1])
    print('maximum-cut projective pairs', len(zvals), 'z-values', sorted(set(zvals.values())))
    print('10.442 lhs/target/failure', lhs, target, target - lhs)
    print('doubled favorable/minimum/failure', favorable_parent_energy,
          child_min_energy, favorable_parent_energy - child_min_energy)


if __name__ == '__main__':
    main()
