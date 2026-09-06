#!/usr/bin/env python3
"""Search failures of max-over-orientation one-step coverage."""

import random
from itertools import product


def energy(A, ids, x):
    return sum(A[i][j] * x[r] * x[s]
               for r, i in enumerate(ids) for s, j in enumerate(ids) if r != s)


def data(A, ids):
    spins = list(product((-1, 1), repeat=len(ids)))
    vals = [energy(A, ids, x) for x in spins]
    P, lo = max(vals), min(vals)
    return P, -lo, [x for x, v in zip(spins, vals) if v == P], [x for x, v in zip(spins, vals) if v == lo]


def pair_coverage(A, p, n):
    size = len(A)
    P, N, _, _ = data(A, tuple(range(size)))
    R = P + N
    total = 0
    rows = []
    for bit in (-1, 1):
        side = tuple(i for i in range(size) if p[i] * n[i] == bit)
        other = tuple(i for i in range(size) if i not in side)
        PX, NX, _, _ = data(A, side)
        QX, IX = max(PX, NX), abs(PX - NX)
        h = energy(A, side, tuple(p[i] for i in side))
        L = sum(abs(sum(A[i][j] * p[i] for i in side)) for j in other)
        mus = (2 * L - (QX - h), 2 * L - (QX + h))
        M = max(0, *mus)
        total += M + IX
        rows.append((side, PX, NX, h, L, mus, M, IX))
    return total, R, rows


random.seed(71059)
for nverts in range(3, 12):
    for trial in range(10000):
        A = [[0] * nverts for _ in range(nverts)]
        for i in range(nverts):
            for j in range(i + 1, nverts):
                A[i][j] = A[j][i] = random.choice((-1, 1))
        P, N, ps, ns = data(A, tuple(range(nverts)))
        pair_results = [(pair_coverage(A, p, n), p, n) for p in ps for n in ns]
        # Report only when every endpoint-pair choice fails the proposed existence claim.
        if max(item[0][0] for item in pair_results) < P + N:
            best = max(pair_results, key=lambda item: item[0][0])
            print("FOUND", nverts, trial, "P,N", P, N, "best", best[0])
            print("p", best[1], "n", best[2])
            for row in A:
                print(" ".join(map(str, row)))
            raise SystemExit
print("none")
