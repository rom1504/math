#!/usr/bin/env python3
"""Finite exact/numerical audit of shared-priority complement consensus."""

from __future__ import annotations

import itertools
import math
import sys
from fractions import Fraction

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from anchored_conflict_r40_check import A5
from envelope_block_cover_r27 import A6, A8, A9


def spins(n: int):
    for tail in itertools.product((-1, 1), repeat=n - 1):
        yield np.asarray((1,) + tail, dtype=np.int64)


def qvalue(a: np.ndarray) -> int:
    return max(abs(int(x @ a @ x)) for x in spins(len(a)))


def complement(a: np.ndarray, s0: tuple[int, ...]) -> np.ndarray:
    b = -a.copy()
    idx = np.asarray(s0)
    b[np.ix_(idx, idx)] = a[np.ix_(idx, idx)]
    return b


def build(a: np.ndarray, m: int, exact_ground: bool, row_cap: int | None = None):
    n = len(a)
    q = qvalue(a)
    selectors = [s for s in itertools.combinations(range(n), m) if 0 in s]
    states = []
    for x in spins(n):
        row = int((a @ x) @ (a @ x))
        for sigma in (-1, 1):
            states.append((sigma, x.copy(), row))
    fibers = {}
    for s0 in selectors:
        b = complement(a, s0)
        values = [sig * int(x @ b @ x) for sig, x, _ in states]
        threshold = max(values) if exact_ground else q
        fibers[s0] = frozenset(i for i, value in enumerate(values)
                                 if value == threshold if exact_ground)
        if not exact_ground:
            fibers[s0] = frozenset(i for i, value in enumerate(values)
                                     if value >= threshold)
        if row_cap is not None:
            fibers[s0] = frozenset(i for i in fibers[s0]
                                     if states[i][2] <= row_cap)
        assert fibers[s0]
    return q, selectors, states, fibers


def minhash_pair(i0, j0, states, coords):
    """Exact expected disagreement and Jaccard-only bound for iid priorities."""
    aset, bset = i0 - j0, j0 - i0
    union = i0 | j0
    total = Fraction(0)
    for i in coords:
        first = int(sum(bool(states[d][1][i] != states[e][1][i])
                        for d in aset for e in j0))
        second = int(sum(bool(states[d][1][i] != states[e][1][i])
                         for d in i0 for e in bset))
        total += Fraction(first, len(union) * len(j0))
        total += Fraction(second, len(union) * len(i0))
    crude = Fraction(len(coords) * (len(union) - len(i0 & j0)), len(union))
    return total, crude


def audit(a: np.ndarray, m: int, name: str, exact_ground: bool,
          row_cap: int | None = None):
    n = len(a)
    q, selectors, states, fibers = build(a, m, exact_ground, row_cap)
    count = len(selectors)
    p = Fraction(m - 1, n - 1)

    # Marginal information and row for the local uniform-fiber kernel.
    marginal = [Fraction(0) for _ in states]
    mean_cond_entropy = 0.0
    row = Fraction(0)
    for s0 in selectors:
        f = fibers[s0]
        mean_cond_entropy += math.log(len(f)) / count
        for d in f:
            marginal[d] += Fraction(1, count * len(f))
            row += Fraction(states[d][2], count * len(f))
    output_entropy = -sum(float(z) * math.log(float(z)) for z in marginal if z)
    information = output_entropy - mean_cond_entropy

    # Independent-kernel conflict from conditional coordinate means.
    independent = Fraction(0)
    for i in range(1, n):
        plus = minus = Fraction(0)
        for s0 in selectors:
            if i not in s0:
                continue
            f = fibers[s0]
            plus += Fraction(int(sum(bool(states[d][1][i] == 1) for d in f)),
                             count * len(f))
            minus += Fraction(int(sum(bool(states[d][1][i] == -1) for d in f)),
                              count * len(f))
        assert plus + minus == p
        independent += 2 * plus * minus / p

    # Shared-priority deterministic-map conflict, exactly averaged over clocks.
    shared = Fraction(0)
    jaccard = Fraction(0)
    pair_colored = []
    for s0 in selectors:
        for t0 in selectors:
            coords = (set(s0) & set(t0)) - {0}
            value, crude = minhash_pair(fibers[s0], fibers[t0], states, coords)
            shared += value / (count * count * p)
            jaccard += crude / (count * count * p)
            if s0 != t0:
                pair_colored.append((value, crude))

    assert shared <= jaccard
    result = {
        "name": name,
        "fiber": "ground" if exact_ground else "incidence",
        "row_cap": row_cap,
        "q": q,
        "selectors": count,
        "fiber_range": (min(map(len, fibers.values())), max(map(len, fibers.values()))),
        "mean_row": float(row),
        "mean_row_exact": str(row),
        "information_nats": information,
        "independent_conflict": float(independent),
        "independent_conflict_exact": str(independent),
        "minhash_conflict": float(shared),
        "minhash_conflict_exact": str(shared),
        "jaccard_bound": float(jaccard),
        "jaccard_bound_exact": str(jaccard),
        "minhash_over_independent": float(shared / independent) if independent else 0.0,
        "colored_over_jaccard": float(shared / jaccard) if jaccard else 0.0,
    }
    print(result)
    return result


def synthetic_formula_check():
    labels = [
        (-1, np.asarray([1, 1, 1]), 0),
        (-1, np.asarray([1, 1, -1]), 0),
        (-1, np.asarray([1, -1, 1]), 0),
        (-1, np.asarray([1, -1, -1]), 0),
    ]
    i0, j0 = frozenset((0, 1, 2)), frozenset((1, 2, 3))
    exact, _ = minhash_pair(i0, j0, labels, {1, 2})
    # Exhaust all priority permutations independently.
    brute = Fraction(0)
    for order in itertools.permutations(range(4)):
        rank = {d: r for r, d in enumerate(order)}
        di, dj = min(i0, key=rank.get), min(j0, key=rank.get)
        brute += Fraction(sum(labels[di][1][i] != labels[dj][1][i] for i in (1, 2)), 24)
    assert exact == brute


def main():
    synthetic_formula_check()
    results = []
    for name, a, m in (("A5", A5, 4), ("A6", A6, 5),
                       ("A8", A8, 6), ("A9", A9, 7)):
        for exact in (False, True):
            results.append(audit(np.asarray(a), m, name, exact))
    for name, a, m, cap in (("A5", A5, 4, 24), ("A6", A6, 5, 30),
                            ("A8", A8, 6, 64), ("A9", A9, 7, 96)):
        results.append(audit(np.asarray(a), m, name, False, cap))
    assert all(r["minhash_conflict"] <= r["jaccard_bound"] + 1e-12 for r in results)
    print("PASS complement_consensus_r42_check")


if __name__ == "__main__":
    main()
