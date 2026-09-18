#!/usr/bin/env python3
"""Audit the sufficient maximum-cut inequality (10.442) on sign matrices."""

import argparse
from itertools import product
from random import Random

from kmin_carleson_r8 import EndpointSystem
from search_kmin_r8 import gauge_matrix


def cuts(A, ids):
    pats = [(1,) + p for p in product((1, -1), repeat=max(0, len(ids)-1))]
    vals = []
    for x in pats:
        val = sum(A[i][j] for u, i in enumerate(ids)
                  for v, j in enumerate(ids) if u < v and x[u] != x[v])
        vals.append(val)
    return pats, vals


def audit(A):
    n = len(A)
    sys = EndpointSystem(A)
    U = (1 << n) - 1
    P, N, ps, ns = sys.endpoints(U)
    worst = None
    checked = 0
    for p in ps:
        for neg in ns:
            signs = tuple(p[i] * neg[i] for i in range(n))
            S = [i for i, q in enumerate(signs) if q == 1]
            T = [i for i, q in enumerate(signs) if q == -1]
            if not S or not T:
                continue
            # Gauge the positive endpoint to all ones.
            G = [[A[i][j] * p[i] * p[j] for j in range(n)] for i in range(n)]
            xp, xc = cuts(G, S)
            yp, yc = cuts(G, T)
            aS = sum(G[i][j] for u, i in enumerate(S) for j in S[u+1:])
            aT = sum(G[i][j] for u, i in enumerate(T) for j in T[u+1:])
            mS, MS = min(xc), max(xc)
            mT, MT = min(yc), max(yc)
            # Only the hard same-positive branch of 10.442.
            if aS + aT < 0 or MS > aS - mS or MT > aT - mT:
                continue
            X = [x for x, q in zip(xp, xc) if q == MS]
            Y = [y for y, q in zip(yp, yc) if q == MT]
            maxz = max(abs(sum(G[i][j] * x[u] * y[v]
                               for u, i in enumerate(S)
                               for v, j in enumerate(T)))
                       for x in X for y in Y)
            target = aS + aT - mS - mT
            slack = 2 * (MS + MT) + maxz - target
            checked += 1
            row = (slack, target, MS + MT, maxz, tuple(S), p, neg,
                   (aS, aT, mS, mT, MS, MT), len(X), len(Y))
            if worst is None or row[0] < worst[0]:
                worst = row
            if slack < 0:
                return checked, row
    return checked, worst


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--n', type=int, default=7)
    ap.add_argument('--trials', type=int, default=0)
    ap.add_argument('--seed', type=int, default=90210)
    args = ap.parse_args()
    total = 1 << ((args.n - 1) * (args.n - 2) // 2)
    rng = Random(args.seed)
    codes = range(total) if not args.trials else [rng.randrange(total) for _ in range(args.trials)]
    count = 0
    global_worst = None
    for z, code in enumerate(codes, 1):
        A = gauge_matrix(args.n, code)
        checked, row = audit(A)
        count += checked
        if row is not None and (global_worst is None or row[0] < global_worst[0]):
            global_worst = (*row, code)
        if row is not None and row[0] < 0:
            print('FAIL code', code, row)
            for r in A:
                print(r)
            return 1
    print('PASS matrices/cases', z, count)
    print('worst', global_worst)


if __name__ == '__main__':
    raise SystemExit(main())
