#!/usr/bin/env python3
"""Exhaustive checks of common nested rearrangement synchronization."""

from __future__ import annotations

from itertools import product
import json


def spin_vectors(n):
    return list(product((-1, 1), repeat=n))


def weight(x):
    return sum(v == 1 for v in x)


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def micro_opt(n, m, h, alpha, beta):
    spins = spin_vectors(n)
    best = None
    for xs in product(spins, repeat=m):
        ks = [weight(x) for x in xs]
        val = sum(h[a][ks[a]] for a in range(m))
        for a in range(m):
            for b in range(a + 1, m):
                sa, sb = 2 * ks[a] - n, 2 * ks[b] - n
                val += alpha[a][b] * dot(xs[a], xs[b])
                val += beta[a][b] * sa * sb
        best = val if best is None else max(best, val)
    return best


def quotient_opt(n, m, h, alpha, beta):
    best = None
    for ks in product(range(n + 1), repeat=m):
        val = sum(h[a][ks[a]] for a in range(m))
        for a in range(m):
            for b in range(a + 1, m):
                sa, sb = 2 * ks[a] - n, 2 * ks[b] - n
                val += alpha[a][b] * (n - 2 * abs(ks[a] - ks[b]))
                val += beta[a][b] * sa * sb
        best = val if best is None else max(best, val)
    return best


def verify():
    checks = 0
    rows = []
    for n, m in ((2, 2), (3, 2), (2, 3), (3, 3), (4, 2)):
        # Deterministic heterogeneous tables and nonnegative identity-channel
        # coefficients; beta has both signs.
        h = [[((a + 2) * k * k + 3 * a * k + n) % 11 - 5
              for k in range(n + 1)] for a in range(m)]
        alpha = [[0 for _ in range(m)] for _ in range(m)]
        beta = [[0 for _ in range(m)] for _ in range(m)]
        for a in range(m):
            for b in range(a + 1, m):
                alpha[a][b] = 1 + (a + 2 * b) % 3
                beta[a][b] = -2 + (2 * a + b) % 5
        microscopic = micro_opt(n, m, h, alpha, beta)
        quotient = quotient_opt(n, m, h, alpha, beta)
        assert microscopic == quotient
        checks += 1
        rows.append({"n": n, "blocks": m, "optimum": microscopic})
    frustration_checks = 0
    for n in (2, 4, 6):
        balanced = [x for x in spin_vectors(n) if weight(x) == n // 2]
        actual = max(dot(x, y) + dot(y, z) - dot(x, z)
                     for x in balanced for y in balanced for z in balanced)
        assert actual == n
        assert 3 * n - actual == 2 * n
        frustration_checks += 1
    unbalanced_cycle_checks = 0
    for ell in (3, 4, 5):
        # One negative edge makes the cycle sign product negative.  Enumerate
        # one-coordinate assignments; the balanced n-copy construction then
        # pairs each maximizer with its global negative.
        signs = [1] * (ell - 1) + [-1]
        one_site = max(sum(signs[a] * z[a] * z[(a + 1) % ell]
                           for a in range(ell))
                       for z in spin_vectors(ell))
        assert one_site == ell - 2
        unbalanced_cycle_checks += 1
    return {"exact_micro_vs_quotient_checks": checks,
            "mixed_sign_frustration_checks": frustration_checks,
            "unbalanced_cycle_checks": unbalanced_cycle_checks,
            "instances": rows}


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
