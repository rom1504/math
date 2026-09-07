#!/usr/bin/env python3
"""Exact-size MILP for a +-1 counterexample to the hard endpoint cut lemma."""

from itertools import combinations, product

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp

from weighted_cut_stopping_lp_r9 import model


def solve(s, t):
    edges, states, C, endpoint, shores = model(s, t)
    # Up to complements and independent relabelling, only the cardinality of
    # each designated child minimum cut matters.
    cases = []
    for rs in range(s // 2 + 1):
        for rt in range(t // 2 + 1):
            pS = (0,) * (s - rs) + (1,) * rs
            pT = (0,) * (t - rt) + (1,) * rt
            i = shores[0][2].index(pS)
            j = shores[1][2].index(pT)
            cases.append((i, j))

    best = None
    ones = np.ones(len(edges))
    for i, j in cases:
        mS, mT = shores[0][1][i], shores[1][1][j]
        rows = []
        rhs = []

        def add(row, bound=0.0):
            # row @ w <= bound with w=2*x-1.
            rows.append(2 * row)
            rhs.append(bound + row @ ones)

        for cut in C:
            add(-cut)                 # every full cut is >= 0
            add(cut - endpoint)       # every full cut is <= c
        for (aX, cuts, _), mX in zip(shores, (mS, mT)):
            for kX in cuts:
                add(mX - kX)          # designated child cut is minimum
                add(kX - aX + mX)     # child is positive-dominant
        a = shores[0][0] + shores[1][0]
        add(-a)                        # a_S+a_T >= 0

        # Maximize a-mS-mT-c, equivalently minimize its negative.
        objective_w = -(a - mS - mT - endpoint)
        objective_x = 2 * objective_w
        res = milp(
            objective_x,
            integrality=np.ones(len(edges)),
            bounds=Bounds(np.zeros(len(edges)), np.ones(len(edges))),
            constraints=LinearConstraint(np.array(rows), -np.inf, np.array(rhs)),
            options={"time_limit": 300, "mip_rel_gap": 0},
        )
        if res.x is None:
            print("case", shores[0][2][i], shores[1][2][j], res.message)
            continue
        w = 2 * np.rint(res.x).astype(int) - 1
        value = int(round((a - mS - mT - endpoint) @ w))
        print("case", shores[0][2][i], shores[1][2][j], "gap", value,
              "status", res.status, "bound", res.mip_gap)
        if best is None or value > best[0]:
            best = (value, w, i, j, edges, shores, C, endpoint)
    return best


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("s", type=int)
    ap.add_argument("t", type=int)
    args = ap.parse_args()
    best = solve(args.s, args.t)
    if best:
        value, w, i, j, edges, shores, C, endpoint = best
        print("BEST GAP", value)
        print("patterns", shores[0][2][i], shores[1][2][j])
        print("weights", list(zip(edges, map(int, w))))
