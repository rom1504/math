#!/usr/bin/env python3
"""Finite checks for the repaired local-field affine-shell statements.

The important repair is the orientation rho: an absolute maximizer x need
not satisfy H_A(x)=Q, and x -> -x does not change a quadratic energy.
"""

from __future__ import annotations

import argparse
import itertools
import random


def spins(n: int):
    return itertools.product((-1, 1), repeat=n)


def energy(edges, x) -> int:
    return sum(a * x[i] * x[j] for i, j, a in edges)


def flip(x, mask: int, indices) -> tuple[int, ...]:
    y = list(x)
    for bit, i in enumerate(indices):
        if mask >> bit & 1:
            y[i] *= -1
    return tuple(y)


def verify(edges, n: int) -> None:
    landscape = [(tuple(x), energy(edges, x)) for x in spins(n)]
    x, hx = max(landscape, key=lambda item: abs(item[1]))
    q = abs(hx)
    rho = 1 if hx >= 0 else -1

    ell = [
        rho * x[i] * sum(a * x[j] for u, j, a in edges if u == i)
        + rho * x[i] * sum(a * x[u] for u, j, a in edges if j == i)
        for i in range(n)
    ]
    assert all(0 <= value <= q for value in ell)
    assert sum(ell) == 2 * q

    for k in range(1, n):
        indices = sorted(range(n), key=lambda i: ell[i])[:k]
        delta = 4 * k * q / n + 2 * k * (k - 1)
        coset = [flip(x, mask, indices) for mask in range(1 << k)]

        assert len(set(coset)) == 1 << k
        representatives = {min(y, tuple(-v for v in y)) for y in coset}
        assert len(representatives) == 1 << k

        for mask, y in enumerate(coset):
            internal = sum(
                rho * a * x[i] * x[j]
                for i, j, a in edges
                if i in indices
                and j in indices
                and (mask >> indices.index(i) & 1)
                and (mask >> indices.index(j) & 1)
            )
            selected_sum = sum(
                ell[i] for bit, i in enumerate(indices) if mask >> bit & 1
            )
            assert rho * energy(edges, y) == q - 2 * selected_sum + 4 * internal
            assert q - abs(energy(edges, y)) <= delta + 1e-12

        # Odd products are exactly XOR in the affine coordinate mask.
        for a in range(1 << k):
            for b in range(1 << k):
                c = (a * 1103515245 + b * 12345) & ((1 << k) - 1)
                product = tuple(
                    coset[a][i] * coset[b][i] * coset[c][i] for i in range(n)
                )
                assert product == coset[a ^ b ^ c]

        if k % 2:
            continue

        ports = [x] + [flip(x, 1 << bit, indices) for bit in range(k)]
        seen_masks = set()
        for epsilon in spins(k + 1):
            field = [sum(e * w[i] for e, w in zip(epsilon, ports)) for i in range(n)]
            assert all(value != 0 for value in field)
            selector = tuple(1 if value > 0 else -1 for value in field)
            t = sum(epsilon)
            oriented = tuple((1 if t > 0 else -1) * value for value in selector)
            selector_mask = 0
            for bit, i in enumerate(indices):
                if oriented[i] != x[i]:
                    selector_mask |= 1 << bit
            assert oriented == coset[selector_mask]
            seen_masks.add(selector_mask)

            # m=0 is exactly B_A(0)=Q; m=1 already checks every endpoint
            # against a full Boolean response maximization.
            for m in (1,):
                g = [m * value for value in field]
                response = max(abs(hy) + sum(gi * y[i] for i, gi in enumerate(g))
                               for y, hy in landscape)
                gap = q + sum(abs(gi) for gi in g) - response
                assert -1e-12 <= gap <= delta + 1e-12

        expected = {0}
        expected.update(
            mask
            for mask in range(1 << k)
            if bin(mask).count("1") in (k // 2, k // 2 + 1)
        )
        assert seen_masks == expected


def edge_list(n: int, signs) -> list[tuple[int, int, int]]:
    pairs = list(itertools.combinations(range(n), 2))
    return [(i, j, a) for (i, j), a in zip(pairs, signs)]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=20260817)
    parser.add_argument("--random-per-order", type=int, default=4)
    args = parser.parse_args()
    rng = random.Random(args.seed)

    checked = 0
    # Exhaust all signings through n=5 (global negations deliberately kept).
    for n in range(2, 6):
        e = n * (n - 1) // 2
        for signs in spins(e):
            verify(edge_list(n, signs), n)
            checked += 1

    # Stress larger orders by random signings, including both global signs.
    for n in range(6, 11):
        e = n * (n - 1) // 2
        for _ in range(args.random_per_order):
            signs = tuple(rng.choice((-1, 1)) for _ in range(e))
            verify(edge_list(n, signs), n)
            verify(edge_list(n, tuple(-a for a in signs)), n)
            checked += 2

    print(f"PASS: {checked} signings; all subsets/endpoints; n <= 10")


if __name__ == "__main__":
    main()
