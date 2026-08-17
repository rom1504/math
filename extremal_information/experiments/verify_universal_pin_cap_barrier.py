#!/usr/bin/env python3
"""Finite exact checks for the universal coordinate-pin cap barrier.

This verifies the adversarial child margin exhaustively for small orders and
checks the rank-one exact compiler's effective landscape.  It is a finite
arithmetic verifier, not evidence for an asymptotic extension beyond UP.1.
"""

from __future__ import annotations

import itertools


def spins(n: int):
    return itertools.product((-1, 1), repeat=n)


def energy(word, x, edges):
    return sum(a * x[i] * x[j] for a, (i, j) in zip(word, edges))


def check_adversarial_margin() -> None:
    for k in range(2, 8):
        u = (1,) * k
        edges = list(itertools.combinations(range(k), 2))
        for x in spins(k):
            d_raw = sum(a != b for a, b in zip(x, u))
            if d_raw > k // 2:
                x = tuple(-z for z in x)
                d_raw = k - d_raw
            S = {i for i in range(k) if x[i] != u[i]}
            word = tuple(-u[i] * u[j] if ((i in S) != (j in S)) else 1
                         for i, j in edges)
            gap = energy(word, x, edges) - energy(word, u, edges)
            assert gap == 2 * d_raw * (k - d_raw)


def rank_one_effective(k: int, u, x) -> int:
    # max_y (x.u)(sum y)+( (sum y)^2-k )/2
    a = sum(x[i] * u[i] for i in range(k))
    return max(a * sum(y) + (sum(y) ** 2 - k) // 2 for y in spins(k))


def check_rank_one_compiler_margin() -> None:
    for k in range(2, 9):
        u = (1,) * k
        gu = rank_one_effective(k, u, u)
        values = []
        for x in spins(k):
            d = min(sum(a != b for a, b in zip(x, u)),
                    sum(a != -b for a, b in zip(x, u)))
            gx = rank_one_effective(k, u, x)
            assert gu - gx == 2 * k * d
            assert gu - gx >= 2 * d * (k - d)
            values.append(gx)
        assert max(values) - min(values) >= 2 * (k // 2) * (k - k // 2)


if __name__ == "__main__":
    check_adversarial_margin()
    check_rank_one_compiler_margin()
    print("universal coordinate-pin cap checks: PASS")
