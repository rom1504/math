#!/usr/bin/env python3
"""MILP search for an actual +/-1 same-positive cut-lemma violation."""

from itertools import combinations, product

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp


def solve(s=5, t=4, min_s=(0, 0, 0, 1, 1), min_t=(0, 0, 0, 0), time=300):
    n = s + t
    edges = list(combinations(range(n), 2))
    idx = {e: k for k, e in enumerate(edges)}

    def row_cut(z, ids=range(n)):
        ids = tuple(ids)
        row = np.zeros(len(edges), dtype=int)
        for i, j in combinations(ids, 2):
            if z[ids.index(i)] != z[ids.index(j)]:
                row[idx[i, j]] = 1
        return row

    endpoint = np.array([int((i < s) != (j < s)) for i, j in edges])
    totals = []
    shore_cuts = []
    mins = []
    for ids, pat in ((tuple(range(s)), min_s), (tuple(range(s, n)), min_t)):
        total = np.array([int(i in ids and j in ids) for i, j in edges])
        cuts = [row_cut(z, ids) for z in product((0, 1), repeat=len(ids))]
        totals.append(total)
        shore_cuts.append(cuts)
        mins.append(row_cut(pat, ids))

    # Work first in edge-sign variables x.  Convert every l <= r.x <= u to
    # binary b via x=2b-1.
    rows, los, his = [], [], []

    def add(r, lo=-np.inf, hi=np.inf):
        rows.append(np.asarray(r, dtype=int))
        const = -int(np.sum(r))
        los.append(lo - const if np.isfinite(lo) else lo)
        his.append(hi - const if np.isfinite(hi) else hi)

    # 0 <= every global cut <= endpoint cut c.
    for z in ((0,) + tail for tail in product((0, 1), repeat=n - 1)):
        cut = row_cut(z)
        add(cut, lo=0)
        add(cut - endpoint, hi=0)

    # Designated child cuts are minima; positive orientation dominates.
    for total, cuts, m in zip(totals, shore_cuts, mins):
        for k in cuts:
            add(m - k, hi=0)
            add(k - total + m, hi=0)

    # a_S+a_T >= 0.
    a = totals[0] + totals[1]
    add(-a, hi=0)

    A = 2 * np.array(rows, dtype=float)
    cons = LinearConstraint(A, np.array(los), np.array(his))
    # Maximize a-mS-mT-c in x variables. Constant is irrelevant to optimizer.
    objx = a - mins[0] - mins[1] - endpoint
    res = milp(
        c=-2 * objx.astype(float),
        integrality=np.ones(len(edges)),
        bounds=Bounds(np.zeros(len(edges)), np.ones(len(edges))),
        constraints=cons,
        options={"time_limit": time, "mip_rel_gap": 0.0},
    )
    print(res.message, "status", res.status, "bound", getattr(res, "mip_dual_bound", None))
    if res.x is not None:
        x = np.rint(2 * res.x - 1).astype(int)
        print("objective", int(objx @ x), "c", int(endpoint @ x),
              "a,mS,mT", int(a @ x), int(mins[0] @ x), int(mins[1] @ x))
        for i in range(n):
            print([0 if i == j else int(x[idx[min(i,j),max(i,j)]]) for j in range(n)])
    return res


if __name__ == "__main__":
    solve()
