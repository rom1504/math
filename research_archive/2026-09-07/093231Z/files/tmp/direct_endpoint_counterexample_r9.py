#!/usr/bin/env python3
"""Targeted exact MILP for an honest signing violating endpoint stopping.

We prescribe p=1, the parent negative endpoint S|T, a trivial positive
child ground on S, and a B|C positive child ground on T.  The MILP maximizes

    h_S + E_T(B|C) + H - 4c,

under 0 <= every parent cut <= c and h_S,h_T,E_T,H >= 0.  A positive
integer optimum is simultaneously a violation of (10.437) and (10.435).
Scratch only; any returned matrix is exhaustively verified in integer
arithmetic.
"""

from __future__ import annotations

import argparse
from itertools import combinations

import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix, vstack


def vectors(n: int, s: int, k: int):
    edges = list(combinations(range(n), 2))
    S = set(range(s))
    T = set(range(s, n))
    C = set(range(n - k, n))

    def coeff(pred):
        return np.array([int(pred(i, j)) for i, j in edges], dtype=np.int8)

    endpoint = coeff(lambda i, j: (i in S) != (j in S))
    intS = coeff(lambda i, j: i in S and j in S)
    intT = coeff(lambda i, j: i in T and j in T)
    childcut = coeff(lambda i, j: i in T and j in T and ((i in C) != (j in C)))
    # Energies in terms of a_e: hX=2 sum internal a, ET=hT-4 childcut.
    hS = 2 * intS
    hT = 2 * intT
    H = hS + hT
    ET = hT - 4 * childcut
    # objective hS+ET+H-4c.
    objective_a = hS + ET + H - 4 * endpoint
    return edges, endpoint, hS, hT, H, ET, objective_a


def signed_constant(v):
    # For a=2x-1, v.a = 2v.x-sum(v).
    return -int(np.sum(v))


def solve(n: int, s: int, k: int, time_limit: float):
    edges, endpoint, hS, hT, H, ET, objective_a = vectors(n, s, k)
    m = len(edges)
    erow = np.flatnonzero(endpoint)
    rows, cols, data, lo, hi = [], [], [], [], []
    rowno = 0
    # Projective cuts have vertex zero outside.  Add all parent endpoint
    # conditions in a sparse matrix.
    for mask in range(1 << (n - 1)):
        crossing = [q for q, (i, j) in enumerate(edges)
                    if (((mask >> (i - 1)) & 1) if i else 0)
                    != (((mask >> (j - 1)) & 1) if j else 0)]
        # cut >= 0: 2 sum_cross x >= |cross|.
        for q in crossing:
            rows.append(rowno); cols.append(q); data.append(2)
        lo.append(len(crossing)); hi.append(np.inf); rowno += 1
        # cut <= endpoint: 2(sum_cross-sum_endpoint)x <= |cross|-|endpoint|.
        mark = {q: 2 for q in crossing}
        for q in erow:
            mark[q] = mark.get(q, 0) - 2
        for q, val in mark.items():
            if val:
                rows.append(rowno); cols.append(q); data.append(val)
        lo.append(-np.inf); hi.append(len(crossing) - len(erow)); rowno += 1
    M = coo_matrix((data, (rows, cols)), shape=(rowno, m)).tocsr()
    # Sign conditions hS,hT,ET >= 0.  H then follows but include it in audit.
    signs = np.vstack((hS, hT, ET))
    M = vstack((M, 2 * signs), format="csr")
    for v in signs:
        lo.append(-signed_constant(v)); hi.append(np.inf)

    # scipy minimizes.  objective_a.a = 2*objective_a.x + constant.
    ans = milp(
        -2.0 * objective_a.astype(float),
        integrality=np.ones(m), bounds=Bounds(0, 1),
        constraints=LinearConstraint(M, np.array(lo), np.array(hi)),
        options={"time_limit": time_limit, "mip_rel_gap": 0.0},
    )
    if ans.x is None:
        return None, ans
    bits = np.rint(ans.x).astype(np.int8)
    avec = 2 * bits - 1
    value = int(objective_a @ avec)
    return (edges, avec, value), ans


def audit(n: int, s: int, k: int, result):
    edges, avec, value = result
    _, endpoint, hS, hT, H, ET, objective_a = vectors(n, s, k)
    c = int(endpoint @ avec)
    mincut, maxcut = 10**9, -10**9
    argmin = argmax = None
    for mask in range(1 << (n - 1)):
        cut = sum(a for (i, j), a in zip(edges, avec)
                  if (((mask >> (i - 1)) & 1) if i else 0)
                  != (((mask >> (j - 1)) & 1) if j else 0))
        if cut < mincut: mincut, argmin = cut, mask
        if cut > maxcut: maxcut, argmax = cut, mask
    vals = {"c": c, "mincut": mincut, "maxcut": maxcut,
            "hS": int(hS @ avec), "hT": int(hT @ avec),
            "H": int(H @ avec), "ET": int(ET @ avec),
            "objective": int(objective_a @ avec),
            "argmin": argmin, "argmax": argmax}
    assert mincut == 0 and maxcut == c
    assert vals["hS"] >= 0 and vals["hT"] >= 0 and vals["ET"] >= 0
    assert vals["objective"] == value
    A = [[0] * n for _ in range(n)]
    for (i, j), a in zip(edges, avec):
        A[i][j] = A[j][i] = int(a)
    return vals, A


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("n", type=int)
    ap.add_argument("s", type=int)
    ap.add_argument("k", type=int, help="size of C in the T=B+C child cut")
    ap.add_argument("--time", type=float, default=120.0)
    args = ap.parse_args()
    result, ans = solve(args.n, args.s, args.k, args.time)
    print("solver", ans.message, "fun", ans.fun, "bound", ans.mip_gap, flush=True)
    if result is None:
        return 1
    vals, A = audit(args.n, args.s, args.k, result)
    print("exact", vals)
    if vals["objective"] > 0:
        print("COUNTEREXAMPLE")
        for row in A:
            print(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
