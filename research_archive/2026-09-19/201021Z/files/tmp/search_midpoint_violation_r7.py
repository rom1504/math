#!/usr/bin/env python3
"""Search endpoint cuts where both child restrictions miss child midpoints."""

import random
from itertools import product


def energy(A, ids, x):
    return sum(A[i][j] * x[r] * x[s]
               for r, i in enumerate(ids) for s, j in enumerate(ids) if r != s)


def endpoint_data(A, ids):
    spins = list(product((-1, 1), repeat=len(ids)))
    vals = [energy(A, ids, x) for x in spins]
    hi, lo = max(vals), min(vals)
    return hi, -lo, [x for x, v in zip(spins, vals) if v == hi], [x for x, v in zip(spins, vals) if v == lo]


random.seed(71058)
for n in range(4, 11):
    for trial in range(3000):
        A = [[0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                A[i][j] = A[j][i] = random.choice((-1, 1))
        P, N, ps, ns = endpoint_data(A, tuple(range(n)))
        found = False
        for p in ps:
            for q in ns:
                sides = [tuple(i for i in range(n) if p[i] * q[i] == bit) for bit in (-1, 1)]
                if not all(sides):
                    continue
                rows = []
                for side in sides:
                    PX, NX, _, _ = endpoint_data(A, side)
                    h = energy(A, side, tuple(p[i] for i in side))
                    rows.append((side, PX, NX, h, abs(PX - NX), 2 * abs(h)))
                if all(row[5] < row[4] for row in rows):
                    print("FOUND", n, trial, "root", P, N, "rows", rows)
                    print("p", p, "n", q)
                    for row in A:
                        print(" ".join(map(str, row)))
                    found = True
                    break
            if found:
                break
        if found:
            raise SystemExit
print("none")
