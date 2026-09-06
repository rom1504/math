#!/usr/bin/env python3
"""Explore conditional block-response games on the exact A9 minimizer."""

from itertools import combinations, product
from math import comb

import numpy as np
from scipy.optimize import linprog

from global_variation_r33_check import A9, qnorm, response_profiles


def integral_score_table(a: np.ndarray, u: tuple[int, ...]):
    edges, d, z, labels = response_profiles(a, u)
    beta = np.asarray(list(product((-1, 1), repeat=len(edges))), dtype=np.int8)
    score = z.astype(np.int16)[None, :] + 2 * beta.astype(np.int16) @ d.astype(np.int16).T
    return score


def continuous_conditional_value(a: np.ndarray, u: tuple[int, ...], mask: np.ndarray):
    edges, d, z, labels = response_profiles(a, u)
    dd = d[mask]
    zz = z[mask]
    c = np.r_[np.zeros(len(edges)), 1.0]
    aub = np.c_[2.0 * dd, -np.ones(len(dd))]
    sol = linprog(c, A_ub=aub, b_ub=-zz,
                  bounds=[(-1.0, 1.0)] * len(edges) + [(None, None)],
                  method="highs")
    assert sol.success
    return float(sol.fun), sol.x[:-1]


def child_data(a: np.ndarray, u: tuple[int, ...]):
    edges, d, z, labels = response_profiles(a, u)
    actual = np.asarray([a[u[i], u[j]] for i, j in edges], dtype=int)
    energy = 2 * d.astype(int) @ actual
    qb = int(energy.max())
    deficit = qb - energy
    return qb, deficit, np.asarray([s for s, x in labels]), energy, z.astype(int)


def main():
    assert qnorm(A9) == 24
    rows = []
    for b in (3, 4, 5, 6):
        for u in combinations(range(9), b):
            qb, deficit, orientation, energy, z = child_data(A9, u)
            score = integral_score_table(A9, u)
            iall = int(score.max(axis=1).min())
            for D in sorted(set(map(int, deficit))):
                good = deficit <= D
                if np.all(good):
                    continue
                cint_bad, _ = continuous_conditional_value(A9, u, ~good)
                cint_good, _ = continuous_conditional_value(A9, u, good)
                ibad = int(score[:, ~good].max(axis=1).min())
                rows.append((b, u, qb, D, int(good.sum()), iall, ibad,
                             round(cint_bad, 8), round(cint_good, 8),
                             int(np.sum(good & (orientation == 1))),
                             int(np.sum(good & (orientation == -1)))))
    # Print strongest exact walls first: all replacements certified by bad states.
    walls = [row for row in rows if row[6] >= 24]
    walls.sort(key=lambda r: (-r[6], r[0], r[1], r[3]))
    print("wall_count", len(walls), "of", len(rows))
    for row in walls[:80]:
        print(row)
    print("continuous_bad_ge_q", sum(row[7] >= 24 - 1e-7 for row in rows))


if __name__ == "__main__":
    main()
