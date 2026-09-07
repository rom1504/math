#!/usr/bin/env python3
"""Exhaustive tests of candidate soft-degeneracy projection inequalities."""

from __future__ import annotations

from itertools import product
import math

import numpy as np
from scipy.special import logsumexp


def spins(n: int) -> np.ndarray:
    return np.asarray([(1,) + t for t in product((-1, 1), repeat=n - 1)], dtype=np.int8)


SPINS = {n: spins(n) for n in range(1, 7)}


def matrix_from_code(n: int, code: int) -> np.ndarray:
    # Switching-normalize the first row to +1.  The remaining edges are free.
    a = np.zeros((n, n), dtype=np.int8)
    a[0, 1:] = 1
    a[1:, 0] = 1
    bit = 0
    for i in range(1, n):
        for j in range(i + 1, n):
            a[i, j] = a[j, i] = 1 if (code >> bit) & 1 else -1
            bit += 1
    return a


def deficit_data(a: np.ndarray, beta: float) -> tuple[int, float]:
    x = SPINS[len(a)]
    e = np.einsum("bi,ij,bj->b", x, a, x, optimize=True).astype(np.int64)
    e = np.concatenate((e, -e))
    q = int(np.max(e))
    return q, float(logsumexp(-beta * (q - e)))


def main() -> None:
    for beta in (0.05, 0.1, 0.2, 0.5, 1.0, 2.0):
        print("beta", beta)
        for n in range(3, 7):
            m = (n - 1) * (n - 2) // 2
            best = (-math.inf, None)
            worst_min_ratio = (-math.inf, None)
            for code in range(1 << m):
                a = matrix_from_code(n, code)
                q, logd = deficit_data(a, beta)
                child = []
                for i in range(n):
                    ids = [j for j in range(n) if j != i]
                    qc, logdc = deficit_data(a[np.ix_(ids, ids)], beta)
                    child.append((q - qc, logdc))
                # Candidate reverse-LW: product D(C_i) <= D(B)^(n-1).
                gap = sum(v[1] for v in child) - (n - 1) * logd
                if gap > best[0]:
                    best = (gap, (code, q, logd, child))
                # Best unpenalized decrement+entropy score = kappa+log2/beta.
                score = max(d + (logd - logdc) / beta for d, logdc in child)
                ratio_gap = 1.5 * q / n - score
                if ratio_gap > worst_min_ratio[0]:
                    worst_min_ratio = (ratio_gap, (code, q, logd, child, score))
            print(" n", n, "max product gap", best[0], "witness", best[1])
            print("    max [1.5Q/n-score]", worst_min_ratio[0], "witness", worst_min_ratio[1])


if __name__ == "__main__":
    main()
