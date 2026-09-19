#!/usr/bin/env python3
"""Exact-cube LP feasibility test for the radial matching-kernel target.

Temporary diagnostic only.  For each small n and fixed lambda, test whether
q_0=1, q_l>=0 can satisfy:
  * K_q(a)>=0 for every a in {+/-1}^E;
  * l q_{l-1}+C(n-2l,2) q_{l+1} >= lambda E q_l.
"""

from __future__ import annotations

import itertools
import math

import numpy as np
from scipy.optimize import linprog


def matchings_by_size(n: int):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    edge_index = {e: k for k, e in enumerate(edges)}
    levels: list[list[int]] = [[] for _ in range(n // 2 + 1)]

    def rec(available: tuple[int, ...], mask: int, size: int) -> None:
        levels[size].append(mask)
        if len(available) < 2:
            return
        i = available[0]
        # Only add an edge incident to the least available vertex at this
        # recursion node; skipping it ensures each partial matching is unique.
        rec(available[1:], mask, size)
        for pos, j in enumerate(available[1:], start=1):
            rest = available[1:pos] + available[pos + 1 :]
            rec(rest, mask | (1 << edge_index[(i, j)]), size + 1)

    rec(tuple(range(n)), 0, 0)
    # The recursion records an existing matching at multiple later skip nodes.
    levels = [sorted(set(level)) for level in levels]
    return edges, levels


def layer_char_sums(n: int) -> np.ndarray:
    edges, levels = matchings_by_size(n)
    count = 1 << len(edges)
    h = np.zeros((count, len(levels)), dtype=np.int32)
    h[:, 0] = 1
    signings = np.arange(count, dtype=np.uint64)
    for ell, masks in enumerate(levels[1:], start=1):
        total = np.zeros(count, dtype=np.int32)
        for mask in masks:
            # chi_M(a)=(-1)^popcount(a & M)
            parity = np.fromiter(
                (bin(int(a) & mask).count("1") & 1 for a in signings),
                dtype=np.int8,
                count=count,
            )
            total += 1 - 2 * parity
        h[:, ell] = total
    return h


def feasible(n: int, lam: float, h: np.ndarray):
    L = n // 2
    E = n * (n - 1) // 2
    aub = []
    bub = []

    # K_q(a)>=0.
    aub.extend((-h).astype(float))
    bub.extend(np.zeros(h.shape[0]))

    # Recurrence, negated for scipy A_ub x <= b_ub.
    for ell in range(L + 1):
        row = np.zeros(L + 1)
        row[ell] += lam * E
        if ell:
            row[ell - 1] -= ell
        if ell < L:
            row[ell + 1] -= math.comb(n - 2 * ell, 2)
        aub.append(row)
        bub.append(0.0)

    result = linprog(
        np.zeros(L + 1),
        A_ub=np.asarray(aub),
        b_ub=np.asarray(bub),
        A_eq=np.asarray([[1.0] + [0.0] * L]),
        b_eq=np.asarray([1.0]),
        bounds=[(0.0, None)] * (L + 1),
        method="highs",
    )
    return result


def threshold(n: int):
    h = layer_char_sums(n)
    lo = 0.0
    hi = 1.01 / math.sqrt(n - 1)
    best = None
    for _ in range(45):
        mid = (lo + hi) / 2
        result = feasible(n, mid, h)
        if result.success:
            lo = mid
            best = result.x
        else:
            hi = mid
    return lo, best, h


def main() -> None:
    for n in range(3, 7):
        lam, q, h = threshold(n)
        ceiling = 1 / math.sqrt(n - 1)
        min_k = float(np.min(h @ q)) if q is not None else float("nan")
        print(
            f"n={n} lambda={lam:.12g} star={ceiling:.12g} "
            f"ratio={lam/ceiling:.8f} minK={min_k:.4g} q={q}"
        )


if __name__ == "__main__":
    main()
