#!/usr/bin/env python3
"""Exact finite checks for the Hadamard synchronization draft."""

from __future__ import annotations

import itertools
import math
import random


def sylvester(level: int):
    W = [[1]]
    for _ in range(level):
        W = [row + row for row in W] + [row + [-v for v in row] for row in W]
    return W


def mv(W, x):
    return [sum(a * b for a, b in zip(row, x)) for row in W]


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def quad_hollow(A, x):
    return sum(A[i][j] * x[i] * x[j] for i in range(len(x)) for j in range(i + 1, len(x)))


def self_dual_bent(m: int):
    q = 1 << m
    out = []
    for u in range(q):
        for v in range(q):
            parity = bin(u & v).count("1") & 1
            out.append(-1 if parity else 1)
    return out


def all_spins(n: int):
    return itertools.product((-1, 1), repeat=n)


def random_filler(n: int, rng: random.Random):
    A = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            A[i][j] = A[j][i] = rng.choice((-1, 1))
    return A


def cap(A):
    return max(abs(quad_hollow(A, x)) for x in all_spins(len(A))) if A else 0


def check_pair_identity(W, trials: int, rng: random.Random):
    k = len(W)
    root = math.isqrt(k)
    assert root * root == k
    checked = 0
    for _ in range(trials):
        x = [rng.choice((-1, 1)) for _ in range(k)]
        y = [rng.choice((-1, 1)) for _ in range(k)]
        z = mv(W, x)
        cross = dot(z, y)
        lhs = k * root - cross
        # Multiply the squared normalized distance by root to stay integral:
        # (root/2)||y-z/root||^2 = sum(root*y-z)^2/(2root).
        numerator = sum((root * yj - zj) ** 2 for yj, zj in zip(y, z))
        assert numerator == 2 * root * lhs
        l1 = sum(abs(v) for v in z)
        sq = sum((abs(v) - root) ** 2 for v in z)
        assert 2 * root * (k * root - l1) == sq
        checked += 1
    return checked


def check_departure(m: int, rng: random.Random):
    level = 2 * m
    W = sylvester(level)
    x0 = self_dual_bent(m)
    k = len(W)
    q = 1 << m
    assert mv(W, x0) == [q * a for a in x0]
    d = max(1, q // 4)
    assert d < q / 2
    S = set(range(d))
    xs = [-a if i in S else a for i, a in enumerate(x0)]
    z = mv(W, xs)
    assert all((1 if a > 0 else -1) == b for a, b in zip(z, x0))
    assert sum(abs(a) for a in z) == k * q - 2 * d * q
    assert dot(z, x0) == k * q - 2 * d * q

    # Build the coherent cut and arbitrary seeded principal fillers.
    A = [[0] * k for _ in range(k)]
    for i in range(k):
        for j in range(i + 1, k):
            if (i in S) != (j in S):
                val = -x0[i] * x0[j]
            else:
                val = rng.choice((-1, 1))
            A[i][j] = A[j][i] = val
    gain = quad_hollow(A, xs) - quad_hollow(A, x0)
    assert gain == 2 * d * (k - d)
    bridge_loss = dot(mv(W, x0), x0) - dot(z, x0)
    assert bridge_loss == 2 * d * q
    assert gain - bridge_loss == 2 * d * (k - d - q)

    # Exact cap check is feasible through k=16 and is evidence only; the
    # theorem uses the analytic filler bound.
    if k <= 16:
        assert cap(A) <= 3 * k * q
    return k


def main():
    rng = random.Random(20260817)
    checks = 0
    for m in (1, 2, 3):
        W = sylvester(2 * m)
        checks += check_pair_identity(W, 64, rng)
        if m >= 2:
            checks += check_departure(m, rng)
    print(f"Hadamard synchronization checks passed: {checks}")


if __name__ == "__main__":
    main()
