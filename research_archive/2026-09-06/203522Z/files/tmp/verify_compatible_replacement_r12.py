#!/usr/bin/env python3
"""Exact finite certificates for Wave 12 compatible-replacement no-go."""

from collections import Counter
from itertools import permutations, product
import numpy as np


def spins(n):
    return np.array([(1,) + z for z in product((-1, 1), repeat=n-1)], dtype=np.int16)


def energies(a):
    x = spins(len(a))
    return np.sum((x @ a) * x, axis=1)


def q(a):
    return int(np.max(np.abs(energies(a))))


def gauge_signing(n, mask):
    a = np.zeros((n, n), dtype=np.int16)
    a[0, 1:] = a[1:, 0] = 1
    k = 0
    for i in range(1, n):
        for j in range(i+1, n):
            a[i, j] = a[j, i] = 1 if (mask >> k) & 1 else -1
            k += 1
    return a


def gauge_minimizers(n):
    best, out = n*(n-1), []
    for mask in range(1 << ((n-1)*(n-2)//2)):
        a = gauge_signing(n, mask)
        qa = q(a)
        if qa < best:
            best, out = qa, [a]
        elif qa == best:
            out.append(a)
    return best, out


def all_minimizers(n):
    qn, reps = gauge_minimizers(n)
    out = {}
    for a in reps:
        for tail in product((-1, 1), repeat=n-1):
            s = np.array((1,) + tail, dtype=np.int16)
            b = a * np.outer(s, s)
            out[b.tobytes()] = b
    return qn, list(out.values()), len(reps)


def layer(c, a):
    y = spins(len(a))
    exposure = 2*np.sum(np.abs(y @ c.T), axis=1)
    return int(max(0, np.max(exposure - q(a) + np.abs(energies(a)))))


def qdiff(a, b):
    return q(a-b)


def split(a, tail):
    tail = tuple(tail)
    head = tuple(i for i in range(len(a)) if i not in tail)
    return a[np.ix_(head, tail)], a[np.ix_(tail, tail)]


A7 = np.array([
    [0, 1, 1, 1, 1, 1, 1],
    [1, 0, 1, 1, -1, -1, 1],
    [1, 1, 0, 1, -1, 1, -1],
    [1, 1, 1, 0, 1, -1, -1],
    [1, -1, -1, 1, 0, -1, -1],
    [1, -1, 1, -1, -1, 0, -1],
    [1, 1, -1, -1, -1, -1, 0],
], dtype=np.int16)

A6 = np.array([
    [0, 1, 1, 1, 1, 1],
    [1, 0, -1, -1, 1, 1],
    [1, -1, 0, 1, -1, 1],
    [1, -1, 1, 0, 1, -1],
    [1, 1, -1, 1, 0, -1],
    [1, 1, 1, -1, -1, 0],
], dtype=np.int16)

DPRIME = np.array([
    [0, -1, -1, 1, 1],
    [-1, 0, 1, 1, -1],
    [-1, 1, 0, -1, 1],
    [1, 1, -1, 0, -1],
    [1, -1, 1, -1, 0],
], dtype=np.int16)


def orbit(a):
    n = len(a)
    out = {}
    for p in permutations(range(n)):
        p = np.array(p)
        ap = a[np.ix_(p, p)]
        for tail in product((-1, 1), repeat=n-1):
            s = np.array((1,) + tail, dtype=np.int16)
            b = ap * np.outer(s, s)
            out[b.tobytes()] = b
    return list(out.values())


if __name__ == '__main__':
    q5, mins5, g5 = all_minimizers(5)
    q6, mins6, g6 = all_minimizers(6)
    q7, reps7 = gauge_minimizers(7)
    g7 = len(reps7)
    assert (q5, q6, q7) == (8, 10, 18)
    assert (g5, g6, g7) == (12, 12, 3240)
    assert q(A7) == q7

    # Singleton boundary: every minimizing replacement raises L by eight.
    c6, d6 = split(A7, (0, 1, 2, 4, 5, 6))
    vals6 = [layer(c6, g) for g in mins6]
    assert (q(d6), layer(c6, d6), min(vals6)) == (18, 0, 8)

    # Two-row boundary: even Q+L rises under the best minimizing replacement.
    c5, d5 = split(A7, (0, 1, 2, 3, 4))
    vals5 = [layer(c5, g) for g in mins5]
    assert (q(d5), layer(c5, d5), min(vals5)) == (12, 4, 12)
    assert q(d5)+layer(c5, d5) == 16
    assert q5+min(vals5) == 20

    # The order-six/order-five example: uniform Q-distance can be 16 while
    # the cross-aware layer difference is zero.
    d = A6[:5, :5]
    c = A6[5:, :5]
    orb = orbit(DPRIME)
    pairs = Counter((qdiff(d, g), layer(c, g)) for g in orb)
    expected = Counter({(0, 2): 1, (8, 6): 5, (8, 10): 10,
                        (16, 2): 11, (16, 6): 55, (16, 10): 110})
    assert len(orb) == 192 and pairs == expected
    assert sum(v*l for (_, l), v in pairs.items()) == 192*33//4

    print('q5,q6,q7 = 8,10,18; gauge minima = 12,12,3240')
    print('singleton child: Q,L,best-minimizer-L = 18,0,8')
    print('2x5 boundary: original Q+L=16, best minimizer Q+L=20')
    print('order-5 orbit distance/layer pairs:', sorted(pairs.items()))
    print('all exact compatible-replacement certificates: PASS')
