#!/usr/bin/env python3
"""Independent exhaustive audit of the order-nine global hybrid wall."""

from collections import Counter
from itertools import product

import numpy as np


A = np.array([
    [0, 1, -1, 1, 1, 1, 1, -1, -1],
    [1, 0, 1, -1, -1, -1, 1, -1, -1],
    [-1, 1, 0, 1, 1, 1, 1, 1, -1],
    [1, -1, 1, 0, 1, 1, 1, -1, 1],
    [1, -1, 1, 1, 0, 1, -1, 1, -1],
    [1, -1, 1, 1, 1, 0, -1, -1, -1],
    [1, 1, 1, 1, -1, -1, 0, 1, 1],
    [-1, -1, 1, -1, 1, -1, 1, 0, -1],
    [-1, -1, -1, 1, -1, -1, 1, -1, 0],
], dtype=np.int16)


def spins(n):
    return np.array([(1,) + t for t in product((-1, 1), repeat=n - 1)], dtype=np.int16)


def energy_vector(a):
    x = spins(len(a))
    return np.sum((x @ a) * x, axis=1)


def q(a):
    return int(np.max(np.abs(energy_vector(a))))


def matrices(n):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for signs in product((-1, 1), repeat=len(edges)):
        a = np.zeros((n, n), dtype=np.int16)
        for (i, j), value in zip(edges, signs):
            a[i, j] = a[j, i] = value
        yield a


def minimizers(n):
    best = n * (n - 1)
    out = []
    for a in matrices(n):
        value = q(a)
        if value < best:
            best, out = value, [a]
        elif value == best:
            out.append(a)
    return best, out


def main():
    left = (0, 1, 2)
    right = tuple(range(3, 9))
    d3 = A[np.ix_(left, left)]
    d6 = A[np.ix_(right, right)]
    cross = A.copy()
    cross[np.ix_(left, left)] = 0
    cross[np.ix_(right, right)] = 0

    assert q(A) == 24
    assert (q(d3), q(d6), q(cross)) == (6, 14, 24)
    q3, gs3 = minimizers(3)
    q6, gs6 = minimizers(6)
    assert (q3, len(gs3)) == (6, 8)
    assert (q6, len(gs6)) == (10, 384)

    distribution = Counter()
    witness = None
    for g3 in gs3:
        for g6 in gs6:
            hybrid = cross.copy()
            hybrid[np.ix_(left, left)] = g3
            hybrid[np.ix_(right, right)] = g6
            value = q(hybrid)
            distribution[value] += 1
            if witness is None or value < witness[0]:
                witness = (value, g3.copy(), g6.copy(), energy_vector(hybrid))

    assert sum(distribution.values()) == 8 * 384
    assert min(distribution) == 28
    assert witness[0] == 28
    assert max(abs(int(v)) for v in witness[3]) == 28
    print("parent/block/cross norms:", q(A), q(d3), q(d6), q(cross))
    print("exact local minima and labelled minimizer counts:", (q3, len(gs3)), (q6, len(gs6)))
    print("hybrid Q distribution over all 3072 replacement pairs:", sorted(distribution.items()))
    print("minimum hybrid Q:", witness[0])
    print("one attaining energy histogram:", sorted(Counter(map(int, witness[3])).items()))
    print("PASS")


if __name__ == "__main__":
    main()
