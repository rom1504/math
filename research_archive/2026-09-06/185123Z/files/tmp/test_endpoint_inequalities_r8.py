#!/usr/bin/env python3
"""Exact endpoint-pair audit for three proposed local inequalities."""

import argparse
from random import Random

from kmin_carleson_r8 import EndpointSystem
from search_kmin_r8 import gauge_matrix


def pair_rows(S, mask, p, n):
    ids = S.ids(mask)
    pos = {i: k for k, i in enumerate(ids)}
    smask = sum(1 << ids[k] for k in range(len(ids)) if p[k] * n[k] == 1)
    tmask = mask ^ smask
    if not smask or not tmask:
        return None
    out = []
    for X, Y in ((smask, tmask), (tmask, smask)):
        xids, yids = S.ids(X), S.ids(Y)
        px = tuple(p[pos[i]] for i in xids)
        h = S.energy(xids, px)
        P, N, _, _ = S.endpoints(X)
        L = sum(abs(sum(S.A[j][i] * p[pos[i]] for i in xids)) for j in yids)
        out.append((X, P, N, max(P, N), h, L))
    return out


def audit(A):
    S = EndpointSystem(A)
    U = (1 << len(A)) - 1
    P, N, ps, ns = S.endpoints(U)
    R = P + N
    bad = [[], [], []]
    margins = [None, None, None]
    pairs = 0
    for p in ps:
        for n in ns:
            rows = pair_rows(S, U, p, n)
            if not rows:
                continue
            pairs += 1
            H = sum(z[4] for z in rows)
            # (i) exact equality L_X=R/4.
            vals = [4 * z[5] - R for z in rows]
            m1 = -max(abs(v) for v in vals)  # zero iff equality; negative is failure.
            if any(vals):
                bad[0].append((p, n, rows, P, N, vals))
            # (ii) sum(Q_X-|h_X|) <= 2 min(P,N).
            m2 = 2 * min(P, N) - sum(z[3] - abs(z[4]) for z in rows)
            if m2 < 0:
                bad[1].append((p, n, rows, P, N, m2))
            # (iii) Q_S+Q_T <= R-|h_S+h_T|.
            m3 = R - abs(H) - sum(z[3] for z in rows)
            if m3 < 0:
                bad[2].append((p, n, rows, P, N, m3))
            for j, m in enumerate((m1, m2, m3)):
                if margins[j] is None or m < margins[j][0]:
                    margins[j] = (m, p, n, rows, P, N)
    return pairs, bad, margins


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=6)
    ap.add_argument("--trials", type=int, default=0, help="0: exhaustive switching gauge")
    ap.add_argument("--seed", type=int, default=91873)
    args = ap.parse_args()
    total = 1 << ((args.n - 1) * (args.n - 2) // 2)
    rng = Random(args.seed)
    codes = range(total) if not args.trials else [rng.randrange(total) for _ in range(args.trials)]
    all_margins = [None, None, None]
    pair_count = 0
    for z, code in enumerate(codes, 1):
        A = gauge_matrix(args.n, code)
        pairs, bad, margins = audit(A)
        pair_count += pairs
        for j in range(3):
            if bad[j]:
                print("FAIL", j + 1, "n/code", args.n, code, bad[j][0])
                for row in A:
                    print(row)
                return
            if all_margins[j] is None or margins[j][0] < all_margins[j][0]:
                all_margins[j] = (margins[j][0], code, margins[j][1:])
    print("PASS matrices/pairs", z, pair_count)
    for j, m in enumerate(all_margins, 1):
        print("minimum margin", j, m)


if __name__ == "__main__":
    main()
