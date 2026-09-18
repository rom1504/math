#!/usr/bin/env python3
"""Weighted cut-cone relaxation of the hard same-positive stopping lemma.

Normalize the endpoint cross cut to c=1.  For fixed shore sizes and
designated child minimum cuts, maximize a_S+a_T-m_S-m_T subject to every
global cut lying in [0,1] and each designated child state being an absolute
positive endpoint.  A value above one falsifies the weighted relaxation.
Scratch only.
"""

from __future__ import annotations

import argparse
from itertools import combinations, product

import numpy as np
from scipy.optimize import linprog


def model(s: int, t: int):
    n = s + t
    edges = list(combinations(range(n), 2))
    edge_index = {e: k for k, e in enumerate(edges)}
    states = [(0,) + z for z in product((0, 1), repeat=n - 1)]
    global_cuts = np.array(
        [[int(z[i] != z[j]) for i, j in edges] for z in states], dtype=float
    )
    endpoint_cut = np.array(
        [int((i < s) != (j < s)) for i, j in edges], dtype=float
    )
    shores = []
    for ids in (list(range(s)), list(range(s, n))):
        total = np.zeros(len(edges))
        for e in combinations(ids, 2):
            total[edge_index[e]] = 1
        cuts = []
        patterns = []
        for tail in product((0, 1), repeat=len(ids) - 1):
            z = (0,) + tail
            row = np.zeros(len(edges))
            for u, v in combinations(range(len(ids)), 2):
                if z[u] != z[v]:
                    row[edge_index[(ids[u], ids[v])]] = 1
            cuts.append(row)
            patterns.append(z)
        shores.append((total, cuts, patterns))
    return edges, states, global_cuts, endpoint_cut, shores


def solve(s: int, t: int, progress=False):
    edges, states, C, endpoint, shores = model(s, t)
    best = (-np.inf, None)
    total_cases = len(shores[0][1]) * len(shores[1][1])
    done = 0
    for i, mS in enumerate(shores[0][1]):
        for j, mT in enumerate(shores[1][1]):
            rows = [*C, *(-C)]
            rhs = [*np.ones(len(C)), *np.zeros(len(C))]
            names = ([('global_upper', z) for z in states]
                     + [('global_lower', z) for z in states])
            for side, ((aX, cuts, patterns), mX) in enumerate(zip(shores, (mS, mT))):
                for kX, pat in zip(cuts, patterns):
                    # Designated cut is minimum.
                    rows.append(mX - kX)
                    rhs.append(0)
                    names.append(('child_min', side, pat))
                    # Positive endpoint dominates every negative child energy:
                    # k_X <= a_X-m_X.
                    rows.append(kX - aX + mX)
                    rhs.append(0)
                    names.append(('child_positive_dominance', side, pat))
                # The hard branch assumes h_X=2a_X has no prescribed sign,
                # only a_S+a_T>=0.  Impose that below jointly.
            a = shores[0][0] + shores[1][0]
            rows.append(-a)
            rhs.append(0)
            names.append(('a_total_nonnegative',))
            objective = -(a - mS - mT)
            res = linprog(
                objective,
                A_ub=np.array(rows), b_ub=np.array(rhs),
                A_eq=np.array([endpoint]), b_eq=np.array([1.0]),
                bounds=[(None, None)] * len(edges), method='highs',
            )
            if res.success and -res.fun > best[0] + 1e-10:
                best = (-res.fun, (i, j, res, names, edges, shores))
            done += 1
            if progress and done % max(1, total_cases // 10) == 0:
                print('progress', done, '/', total_cases, 'best', best[0], flush=True)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('s', type=int)
    ap.add_argument('t', type=int)
    ap.add_argument('--progress', action='store_true')
    ap.add_argument('--dual', action='store_true')
    args = ap.parse_args()
    value, data = solve(args.s, args.t, args.progress)
    print('sizes', args.s, args.t, 'optimum', value)
    if data is None:
        return 1
    i, j, res, names, edges, shores = data
    print('min_patterns', shores[0][2][i], shores[1][2][j])
    a = shores[0][0] + shores[1][0]
    mS, mT = shores[0][1][i], shores[1][1][j]
    print('a,mS,mT', a @ res.x, mS @ res.x, mT @ res.x)
    print('weights', [(e, round(float(x), 12)) for e, x in zip(edges, res.x)])
    if args.dual:
        print('dual_nonzero')
        for name, dual, slack in zip(names, res.ineqlin.marginals, res.ineqlin.residual):
            if abs(dual) > 1e-9:
                print(name, 'dual', dual, 'slack', slack)
        print('endpoint_dual', res.eqlin.marginals[0])
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
