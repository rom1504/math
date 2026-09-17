#!/usr/bin/env python3
"""Seek a strict 4+5 weighted endpoint counterexample.

Scratch only.  The endpoint cut is normalized to one.  Every other
nontrivial parent cut is required to lie in [delta,1-delta], while fixed
child cuts realize the two positive-dominant child optima.  All output is
independently rechecked after rational reconstruction.
"""

from fractions import Fraction
from itertools import combinations, product

import numpy as np
from scipy.optimize import linprog


def solve(epsilon: Fraction):
    s, t = 4, 5
    n = s + t
    edges = list(combinations(range(n), 2))
    ei = {e: k for k, e in enumerate(edges)}
    states = [(0,) + z for z in product((0, 1), repeat=n - 1)]

    def cut_row(z):
        return np.array([int(z[i] != z[j]) for i, j in edges], dtype=float)

    C = [cut_row(z) for z in states]
    endpoint_state = (0,) * s + (1,) * t
    endpoint = cut_row(endpoint_state)

    shore_ids = (tuple(range(s)), tuple(range(s, n)))
    totals = []
    child_cuts = []
    for ids in shore_ids:
        total = np.zeros(len(edges))
        for e in combinations(ids, 2):
            total[ei[e]] = 1
        cuts = []
        for tail in product((0, 1), repeat=len(ids) - 1):
            z = (0,) + tail
            row = np.zeros(len(edges))
            for u, v in combinations(range(len(ids)), 2):
                if z[u] != z[v]:
                    row[ei[(ids[u], ids[v])]] = 1
            cuts.append(row)
        totals.append(total)
        child_cuts.append(cuts)

    # The non-strict optimum uses the trivial S minimum and the B|C cut in T.
    mS = child_cuts[0][0]
    mT = child_cuts[1][3]  # local pattern 00011
    a = totals[0] + totals[1]

    # Variables are edge weights followed by delta.
    rows, rhs = [], []
    for z, row in zip(states, C):
        if z == (0,) * n or z == endpoint_state:
            continue
        # delta <= cut and cut <= 1-delta.
        rows.append(np.r_[-row, 1.0])
        rhs.append(0.0)
        rows.append(np.r_[row, 1.0])
        rhs.append(1.0)

    for total, cuts, chosen in zip(totals, child_cuts, (mS, mT)):
        for row in cuts:
            # chosen is a child minimum.
            rows.append(np.r_[chosen - row, 0.0])
            rhs.append(0.0)
            # positive dominance: every child cut <= a_X-chosen.
            rows.append(np.r_[row - total + chosen, 0.0])
            rhs.append(0.0)
    rows.append(np.r_[-a, 0.0])
    rhs.append(0.0)
    # a-mS-mT >= 1+epsilon.
    rows.append(np.r_[-a + mS + mT, 0.0])
    rhs.append(-1.0 - float(epsilon))

    objective = np.zeros(len(edges) + 1)
    objective[-1] = -1.0
    ans = linprog(
        objective,
        A_ub=np.array(rows), b_ub=np.array(rhs),
        A_eq=np.array([np.r_[endpoint, 0.0]]), b_eq=np.array([1.0]),
        bounds=[(None, None)] * len(edges) + [(0, None)], method="highs",
    )
    if not ans.success:
        return None
    fr = [Fraction(float(x)).limit_denominator(10**6) for x in ans.x]
    return edges, states, endpoint_state, shore_ids, fr, epsilon


def audit(data):
    edges, states, endpoint_state, shore_ids, vals, epsilon = data
    weights, delta = vals[:-1], vals[-1]

    def cut(z, ids=None):
        allowed = set(ids) if ids is not None else None
        return sum(w for (i, j), w in zip(edges, weights)
                   if (allowed is None or (i in allowed and j in allowed))
                   and z[i] != z[j])

    cuts = [(z, cut(z)) for z in states]
    endpoint_value = dict(cuts)[endpoint_state]
    interior = [v for z, v in cuts if z != (0,) * 9 and z != endpoint_state]
    child = []
    designated = ((0, 0, 0, 0), (0, 0, 0, 1, 1))
    for ids, dz in zip(shore_ids, designated):
        total = sum(w for (i, j), w in zip(edges, weights) if i in ids and j in ids)
        local = []
        for tail in product((0, 1), repeat=len(ids) - 1):
            zz = (0,) + tail
            z = [0] * 9
            for i, b in zip(ids, zz):
                z[i] = b
            local.append((zz, cut(tuple(z), ids)))
        m = dict(local)[dz]
        child.append((total, m, min(v for _, v in local), max(v for _, v in local)))
    objective = sum(x[0] - x[1] for x in child)
    assert endpoint_value == 1
    assert min(interior) >= delta and max(interior) <= 1 - delta
    assert all(m == lo and hi <= total - m for total, m, lo, hi in child)
    assert objective >= 1 + epsilon
    return delta, min(interior), max(interior), objective, child, weights


def main():
    for epsilon in map(Fraction, ("1/48", "1/96", "1/192", "1/1000")):
        data = solve(epsilon)
        if data is None:
            print("epsilon", epsilon, "infeasible")
            continue
        result = audit(data)
        print("epsilon", epsilon, "delta/min/max/objective/children",
              result[:-1])
        print("weights", list(zip(data[0], result[-1])))


if __name__ == "__main__":
    main()
