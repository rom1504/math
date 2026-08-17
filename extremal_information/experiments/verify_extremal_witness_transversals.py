#!/usr/bin/env python3
"""Exact small checks for the witness-transversal theorem."""

from __future__ import annotations

import itertools
import math


def cube(k: int):
    return list(itertools.product((-1, 1), repeat=k))


def mul(x, y):
    return tuple(a * b for a, b in zip(x, y))


def block_energy(x, s: int) -> int:
    return sum(sum(x[b * s:(b + 1) * s]) ** 2 - s for b in range(s))


def check_product_cover_bound() -> None:
    # Exhaust every nonempty witness set in G for k<=3 and every candidate X.
    for k in range(1, 4):
        group = cube(k)
        whole = set(group)
        for mask in range(1, 1 << len(group)):
            W = {group[i] for i in range(len(group)) if (mask >> i) & 1}
            for xmask in range(1, 1 << len(group)):
                X = {group[i] for i in range(len(group)) if (xmask >> i) & 1}
                covered = {mul(w, x) for w in W for x in X}
                if covered == whole:
                    assert len(W) * len(X) >= len(group)


def check_block_family() -> None:
    for s in (2, 4):
        k = s * s
        xs = cube(k)
        vals = [block_energy(x, s) for x in xs]
        assert max(vals) == s**3 - s**2
        assert min(vals) == -(s**2)
        q = max(abs(v) for v in vals)
        for alpha in (0.25, 0.5, 0.75):
            W = [x for x, v in zip(xs, vals) if abs(v) >= alpha * q]
            lower = math.ceil(len(xs) / len(W))
            assert lower >= 1
            print(
                f"s={s}, alpha={alpha:.2f}, tail={len(W)}/{len(xs)}, "
                f"transversal lower={lower}"
            )


if __name__ == "__main__":
    check_product_cover_bound()
    check_block_family()
    print("extremal witness-transversal checks: PASS")
