#!/usr/bin/env python3
"""Exact finite audit of the cross/internal-excess Pareto data for A9.

All quadratic norms are evaluated over projective Boolean spins.  The script
enumerates every set partition of nine vertices, so it also catches the
singleton and whole-block tautologies explicitly.
"""

from collections import defaultdict
from functools import lru_cache
from itertools import product

import numpy as np


A = np.array(
    [
        [0, 1, -1, 1, 1, 1, 1, -1, -1],
        [1, 0, 1, -1, -1, -1, 1, -1, -1],
        [-1, 1, 0, 1, 1, 1, 1, 1, -1],
        [1, -1, 1, 0, 1, 1, 1, -1, 1],
        [1, -1, 1, 1, 0, 1, -1, 1, -1],
        [1, -1, 1, 1, 1, 0, -1, -1, -1],
        [1, 1, 1, 1, -1, -1, 0, 1, 1],
        [-1, -1, 1, -1, 1, -1, 1, 0, -1],
        [-1, -1, -1, 1, -1, -1, 1, -1, 0],
    ],
    dtype=np.int64,
)

QMIN = {1: 0, 2: 2, 3: 6, 4: 8, 5: 8, 6: 10, 7: 18, 8: 20, 9: 24}


def qnorm(M: np.ndarray) -> int:
    n = len(M)
    spins = np.array([(1,) + t for t in product((-1, 1), repeat=n - 1)], dtype=np.int64)
    energies = np.einsum("bi,ij,bj->b", spins, M, spins)
    return int(np.max(np.abs(energies)))


def gauge_signing(n: int, mask: int) -> np.ndarray:
    """First-row-positive representative of a switching class."""
    G = np.zeros((n, n), dtype=np.int64)
    G[0, 1:] = G[1:, 0] = 1
    bit = 0
    for i in range(1, n):
        for j in range(i + 1, n):
            G[i, j] = G[j, i] = 1 if (mask >> bit) & 1 else -1
            bit += 1
    return G


@lru_cache(None)
def all_minimizers(n: int):
    """Enumerate all labelled minimizers, by gauge reps then switchings."""
    best = n * (n - 1)
    reps = []
    for mask in range(1 << ((n - 1) * (n - 2) // 2)):
        G = gauge_signing(n, mask)
        value = qnorm(G)
        if value < best:
            best, reps = value, [G]
        elif value == best:
            reps.append(G)
    out = {}
    for G in reps:
        for tail in product((-1, 1), repeat=n - 1):
            switch = np.array((1,) + tail, dtype=np.int64)
            H = G * np.outer(switch, switch)
            out[H.tobytes()] = H
    return best, tuple(out.values())


def set_partitions(n: int):
    """Yield each set partition in restricted-growth-string order."""
    labels = [0] * n

    def rec(i: int, largest: int):
        if i == n:
            blocks = [[] for _ in range(largest + 1)]
            for v, label in enumerate(labels):
                blocks[label].append(v)
            yield tuple(tuple(block) for block in blocks)
            return
        for label in range(largest + 2):
            labels[i] = label
            yield from rec(i + 1, max(largest, label))

    yield from rec(1, 0)


def data(blocks):
    C = A.copy()
    S = 0
    B = 0
    for block in blocks:
        ix = np.ix_(block, block)
        D = A[ix]
        S += qnorm(D)
        B += QMIN[len(block)]
        C[ix] = 0
    qc = qnorm(C)
    X = S - B
    O = max(qc - QMIN[9], 0)
    # The strongest scalar fact from local minimization plus triangle.
    assert QMIN[9] <= qc + B
    # Common-mosaic response interval for every minimizing replacement.
    assert X <= X + B + qc - QMIN[9]
    return qc, O, X, B, S


def main():
    assert qnorm(A) == QMIN[9]
    rows = []
    profiles = defaultdict(list)
    for blocks in set_partitions(9):
        row = data(blocks)
        sizes = tuple(sorted(map(len, blocks)))
        rows.append((row, sizes, blocks))
        profiles[sizes].append((row, blocks))

    # Non-dominated means no partition has weakly smaller O and weakly larger X.
    frontier = []
    for row, sizes, blocks in rows:
        _, O, X, _, _ = row
        if not any(
            row2[1] <= O and row2[2] >= X and (row2[1] < O or row2[2] > X)
            for row2, _, _ in rows
        ):
            frontier.append((O, X, row[0], row[3], row[4], sizes, blocks))

    print(f"partitions={len(rows)} profiles={len(profiles)}")
    print("global Pareto (O,X,Q(C),B,S,sizes), duplicates suppressed:")
    seen = set()
    for O, X, qc, B, S, sizes, blocks in sorted(frontier):
        key = (O, X, qc, B, S, sizes)
        if key not in seen:
            seen.add(key)
            print(key, blocks)

    print("profile summaries: sizes minO maxX_at_minO maxX")
    for sizes, items in sorted(profiles.items()):
        if sizes in {(1,) * 9, (9,), (3, 6), (1, 8), (2, 7), (4, 5), (3, 3, 3)}:
            min_o = min(row[1] for row, _ in items)
            max_x_at = max(row[2] for row, _ in items if row[1] == min_o)
            max_x = max(row[2] for row, _ in items)
            print(sizes, min_o, max_x_at, max_x)

    a9_seed = next(
        row
        for row, sizes, blocks in rows
        if blocks == ((0, 1, 2), (3, 4, 5, 6, 7, 8))
    )
    assert a9_seed == (24, 0, 4, 16, 20), a9_seed
    print("A9 3+6 seed (Q(C),O,X,B,S)=", a9_seed)

    # Stronger wall: perfect overshoot, large capture, and an exact-minimizer
    # hybrid coexist, but X cannot be subtracted from the scalar triangle bound.
    wall_blocks = ((0, 1, 2, 4, 5, 6), (3,), (7,), (8,))
    wall = data(wall_blocks)
    assert wall == (22, 0, 12, 10, 22), wall
    C = A.copy()
    for block in wall_blocks:
        C[np.ix_(block, block)] = 0
    best = 10**9
    count = 0
    q6, minimizers6 = all_minimizers(6)
    assert q6 == QMIN[6] and len(minimizers6) == 384
    for G in minimizers6:
        hybrid = C.copy()
        hybrid[np.ix_(wall_blocks[0], wall_blocks[0])] = G
        value = qnorm(hybrid)
        if value < best:
            best, count = value, 1
        elif value == best:
            count += 1
    assert best == QMIN[9] and count == 40
    qc, O, X, B, S = wall
    assert qc + B - X == 20 < QMIN[9]
    print("capture wall (Q(C),O,X,B,S,best hybrid,count)=", wall, best, count)


if __name__ == "__main__":
    main()
