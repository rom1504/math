#!/usr/bin/env python3
"""Exact finite audit of block-constant approximate-envelope codebooks.

All quadratic quantities except p^(3/2) q are integer/rational.  We enumerate
oriented projective cuts, selectors, and fixed set partitions exactly; reported
residuals are floating evaluations of the resulting exact expressions.
"""

from __future__ import annotations

import itertools
import math
from collections import Counter

import numpy as np


A6 = np.array([
    [0, 1, 1, 1, 1, 1],
    [1, 0, -1, -1, 1, 1],
    [1, -1, 0, 1, -1, 1],
    [1, -1, 1, 0, 1, -1],
    [1, 1, -1, 1, 0, -1],
    [1, 1, 1, -1, -1, 0],
], dtype=np.int64)

A8 = np.array([
    [0, 1, 1, 1, 1, 1, 1, 1],
    [1, 0, 1, -1, 1, 1, -1, -1],
    [1, 1, 0, 1, -1, 1, -1, -1],
    [1, -1, 1, 0, -1, -1, -1, 1],
    [1, 1, -1, -1, 0, -1, 1, -1],
    [1, 1, 1, -1, -1, 0, 1, 1],
    [1, -1, -1, -1, 1, 1, 0, 1],
    [1, -1, -1, 1, -1, 1, 1, 0],
], dtype=np.int64)

A9 = np.array([
    [0, 1, -1, 1, 1, 1, 1, -1, -1],
    [1, 0, 1, -1, -1, -1, 1, -1, -1],
    [-1, 1, 0, 1, 1, 1, 1, 1, -1],
    [1, -1, 1, 0, 1, 1, 1, -1, 1],
    [1, -1, 1, 1, 0, 1, -1, 1, -1],
    [1, -1, 1, 1, 1, 0, -1, -1, -1],
    [1, 1, 1, 1, -1, -1, 0, 1, 1],
    [-1, -1, 1, -1, 1, -1, 1, 0, -1],
    [-1, -1, -1, 1, -1, -1, 1, -1, 0],
], dtype=np.int64)


def projective_spins(n: int):
    for tail in itertools.product((-1, 1), repeat=n - 1):
        yield np.array((1,) + tail, dtype=np.int64)


def set_partitions_with_sizes(n: int, sizes: tuple[int, ...]):
    """Generate each unlabeled set partition with the prescribed block sizes."""
    # Put exceptional small blocks first.  For each following run of equal
    # sizes, anchoring the least remaining point quotients permutations of
    # those equal blocks without accidentally forcing that point out of a
    # later exceptional block (e.g. the singleton in 2+2+2+2+1).
    sizes = tuple(sorted(sizes))

    def rec(remaining: tuple[int, ...], todo: tuple[int, ...], blocks):
        if not todo:
            yield tuple(tuple(b) for b in blocks)
            return
        size = todo[0]
        # If equal-sized blocks remain, force the least remaining point into
        # this block to quotient permutations of those blocks.
        equal_next = len(todo) > 1 and todo[1] == size
        anchor = remaining[0] if equal_next else None
        pool = remaining[1:] if anchor is not None else remaining
        choose = size - 1 if anchor is not None else size
        for rest in itertools.combinations(pool, choose):
            block = ((anchor,) + rest) if anchor is not None else rest
            block_set = set(block)
            nxt = tuple(i for i in remaining if i not in block_set)
            yield from rec(nxt, todo[1:], blocks + [tuple(sorted(block))])

    yield from rec(tuple(range(n)), sizes, [])


def all_data(A: np.ndarray, m: int):
    n = len(A)
    xs = list(projective_spins(n))
    raw = np.array([int(x @ A @ x) for x in xs], dtype=np.int64)
    q = int(np.max(np.abs(raw)))
    selectors = list(itertools.combinations(range(n), m))
    qS = []
    for S in selectors:
        B = A[np.ix_(S, S)]
        qS.append(max(abs(int(y @ B @ y)) for y in projective_spins(m)))
    qS = np.array(qS, dtype=np.int64)
    p2_num, p2_den = m * (m - 1), n * (n - 1)

    records = []
    for x, e0 in zip(xs, raw):
        row_cost = int((A @ x) @ (A @ x))
        child0 = np.array([
            int(x[list(S)] @ A[np.ix_(S, S)] @ x[list(S)])
            for S in selectors
        ], dtype=np.int64)
        for sigma in (-1, 1):
            # p2_den * X_d(S), exactly.
            xnum = p2_den * sigma * child0 - p2_num * sigma * e0
            records.append((tuple(int(v) for v in x), sigma, row_cost, xnum))
    return q, selectors, records, qS, p2_den


def block_word_keys(partition):
    keys = set()
    k = len(partition)
    # Fix the first block sign to +1 (projective quotient); retain orientation.
    for tail in itertools.product((-1, 1), repeat=k - 1):
        block_signs = (1,) + tail
        x = [0] * sum(len(b) for b in partition)
        for sign, block in zip(block_signs, partition):
            for i in block:
                x[i] = sign
        for sigma in (-1, 1):
            keys.add((tuple(x), sigma))
    assert len(keys) == 2**k
    return keys


def metrics(records, qS, den, q, m, n, allowed_keys, budget):
    eligible = [r for r in records if r[2] <= budget and (r[0], r[1]) in allowed_keys]
    if not eligible:
        return math.inf, len(qS), 0, (), (10**100, den)
    max_xnum = np.stack([r[3] for r in eligible], axis=0).max(axis=0)
    # base_num/den = Q(A[S])-max_d X_d(S), exactly.  Subtracting
    # q(m/n)^(3/2) gives the residual in (10.808).
    base_num = den * qS - max_xnum
    worst_num = int(base_num.max())
    scale = q * (m / n) ** 1.5
    worst_t = worst_num / den - scale
    uncovered = 0
    for a in base_num:
        aa = int(a)
        if aa > 0 and aa * aa * n**3 > q * q * m**3 * den * den:
            uncovered += 1
    return (
        float(worst_t),
        uncovered,
        len(eligible),
        tuple(sorted(Counter(r[2] for r in eligible).items())),
        (worst_num, den),
    )


def audit(A, name, sizes, budgets, ms):
    n = len(A)
    partitions = list(set_partitions_with_sizes(n, sizes))
    print(f"{name} sizes={sizes} partitions={len(partitions)} words={2**len(sizes)}")
    for m in ms:
        q, selectors, records, qS, den = all_data(A, m)
        all_keys = {(r[0], r[1]) for r in records}
        for budget in budgets:
            full = metrics(records, qS, den, q, m, n, all_keys, budget)
            best = None
            best_part = None
            for partition in partitions:
                got = metrics(records, qS, den, q, m, n, block_word_keys(partition), budget)
                # The irrational scale term is common, so exact integer
                # worst-base comparison chooses the optimum partition.
                key = (got[4][0], got[1], -got[2])
                if best is None or key < (best[4][0], best[1], -best[2]):
                    best, best_part = got, partition
            print({
                "m": m,
                "q": q,
                "budget": budget,
                "unrestricted": {"worst_t": full[0], "worst_base": full[4], "uncovered_t0": full[1], "eligible": full[2]},
                "best_block": {"worst_t": best[0], "worst_base": best[4], "uncovered_t0": best[1], "eligible": best[2]},
                "partition": best_part,
                "cost_hist": best[3],
            })


def verify_flip_identities(A):
    n = len(A)
    rng = np.random.default_rng(2702)
    for _ in range(200):
        x = rng.choice((-1, 1), n).astype(np.int64)
        sigma = int(rng.choice((-1, 1)))
        H = set(np.flatnonzero(rng.integers(0, 2, n)).tolist())
        xp = x.copy()
        if H:
            xp[list(H)] *= -1
        S = set(np.flatnonzero(rng.integers(0, 2, n)).tolist())
        shore_full = sum(int(A[i, j] * x[i] * x[j]) for i in H for j in range(n) if j not in H)
        shore_S = sum(int(A[i, j] * x[i] * x[j]) for i in H & S for j in S if j not in H)
        e, ep = sigma * int(x @ A @ x), sigma * int(xp @ A @ xp)
        c = sigma * int(x[list(S)] @ A[np.ix_(list(S), list(S))] @ x[list(S)]) if S else 0
        cp = sigma * int(xp[list(S)] @ A[np.ix_(list(S), list(S))] @ xp[list(S)]) if S else 0
        assert ep - e == -4 * sigma * shore_full
        assert cp - c == -4 * sigma * shore_S
        old_r2, new_r2 = int((A @ x) @ (A @ x)), int((A @ xp) @ (A @ xp))
        AHx = A[:, list(H)] @ x[list(H)] if H else np.zeros(n, dtype=np.int64)
        rhs = -4 * sum(int(x[i] * (A @ A @ x)[i]) for i in H) + 4 * int(AHx @ AHx)
        assert new_r2 - old_r2 == rhs


def main():
    for A in (A6, A8, A9):
        verify_flip_identities(A)
    # Pair/triple partitions are the smallest finite analogues of the
    # entropy-matched block-constant proposal.  Budgets include the old
    # one-deletion audit caps plus n^2 and a looser reference cap.
    audit(A6, "A6", (2, 2, 2), (30, 36, 60), (3, 4, 5))
    audit(A8, "A8", (2, 2, 2, 2), (40, 64, 80), (4, 6, 7))
    audit(A9, "A9-pairs", (2, 2, 2, 2, 1), (40, 80, 112), (5, 7, 8))
    audit(A9, "A9-triples", (3, 3, 3), (40, 80, 112), (5, 7, 8))
    print("PASS: flip identities and exact finite block-cover enumeration")


if __name__ == "__main__":
    main()
