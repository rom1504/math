#!/usr/bin/env python3
"""Random/local search for endpoint-cut inequalities proposed in Wave 8.

Scratch only.  Uses doubled x^T A x normalization and projective spins.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from itertools import product
import random

import numpy as np


def projective_spins(n: int) -> np.ndarray:
    if n == 0:
        return np.zeros((1, 0), dtype=np.int16)
    tails = np.array(list(product((-1, 1), repeat=n - 1)), dtype=np.int16)
    return np.column_stack((np.ones(len(tails), dtype=np.int16), tails))


def extrema(A: np.ndarray, X: np.ndarray):
    vals = np.einsum("bi,ij,bj->b", X, A, X, optimize=True)
    hi, lo = int(vals.max(initial=0)), int(vals.min(initial=0))
    return hi, -lo, X[vals == hi], X[vals == lo]


@dataclass
class Audit:
    q_slack: int
    cap_slack: int
    data: tuple


def audit_matrix(A: np.ndarray, spin_cache: dict[int, np.ndarray]) -> Audit:
    n = len(A)
    X = spin_cache.setdefault(n, projective_spins(n))
    P, N, ps, ns = extrema(A, X)
    R = P + N
    best_q = 10**18
    best_cap = 10**18
    best_data = None
    for p in ps:
        for neg in ns:
            same = p * neg == 1
            if same.all() or (~same).all():
                continue
            pair = []
            cap_sum = 0
            for mask in (same, ~same):
                other = ~mask
                AX = A[np.ix_(mask, mask)]
                XX = spin_cache.setdefault(int(mask.sum()), projective_spins(int(mask.sum())))
                PX, NX, _, _ = extrema(AX, XX)
                QX = max(PX, NX)
                px = p[mask]
                h = int(px @ AX @ px)
                fields = A[np.ix_(other, mask)] @ px
                L = int(np.abs(fields).sum())
                mu_plus = max(0, 2 * L - (QX - h))
                mu_minus = max(0, 2 * L - (QX + h))
                cap_sum += max(mu_plus, mu_minus)
                pair.append((PX, NX, QX, h, L, mu_plus, mu_minus))
            H = pair[0][3] + pair[1][3]
            I = abs(P - N)
            q_slack = R - abs(H) - pair[0][2] - pair[1][2]
            cap_slack = cap_sum - I
            if q_slack < best_q or cap_slack < best_cap:
                best_data = (P, N, tuple(map(tuple, ps)), tuple(map(tuple, ns)),
                             tuple(pair), H, cap_sum, I, A.copy())
            best_q = min(best_q, q_slack)
            best_cap = min(best_cap, cap_slack)
    if best_data is None:
        raise AssertionError("no nontrivial endpoint cut")
    return Audit(best_q, best_cap, best_data)


def random_signing(n: int, rng: random.Random) -> np.ndarray:
    A = np.zeros((n, n), dtype=np.int16)
    for i in range(n):
        for j in range(i + 1, n):
            A[i, j] = A[j, i] = rng.choice((-1, 1))
    return A


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-min", type=int, default=7)
    ap.add_argument("--n-max", type=int, default=13)
    ap.add_argument("--trials", type=int, default=300)
    ap.add_argument("--seed", type=int, default=824691)
    ap.add_argument("--hill-steps", type=int, default=0,
                    help="if positive, run sign-flip simulated annealing instead")
    ap.add_argument("--restarts", type=int, default=4)
    args = ap.parse_args()
    rng = random.Random(args.seed)
    cache = {}
    if args.hill_steps:
        # Lexicographically target the necessary capacity inequality first,
        # then the stronger Q inequality.  A small temperature permits escape
        # from endpoint plateaux while retaining exact integer scoring.
        for n in range(args.n_min, args.n_max + 1):
            best = (10**18, 10**18, None)
            for restart in range(args.restarts):
                A = random_signing(n, rng)
                cur = audit_matrix(A, cache)
                cur_score = (cur.cap_slack, cur.q_slack)
                for step in range(args.hill_steps):
                    i = rng.randrange(n)
                    j = rng.randrange(n - 1)
                    if j >= i:
                        j += 1
                    if i > j:
                        i, j = j, i
                    A[i, j] *= -1
                    A[j, i] *= -1
                    nxt = audit_matrix(A, cache)
                    nxt_score = (nxt.cap_slack, nxt.q_slack)
                    temp = max(0.05, 2.0 * (1.0 - step / args.hill_steps))
                    delta = (nxt_score[0] - cur_score[0]) + 0.05 * (nxt_score[1] - cur_score[1])
                    accept = delta <= 0 or rng.random() < np.exp(-delta / temp)
                    if accept:
                        cur, cur_score = nxt, nxt_score
                    else:
                        A[i, j] *= -1
                        A[j, i] *= -1
                    candidate = (cur.cap_slack, cur.q_slack, cur.data)
                    if candidate[:2] < best[:2]:
                        best = candidate
                        print("HILL", n, restart, step, "cap", best[0], "q", best[1], flush=True)
                    if cur.cap_slack < 0 or cur.q_slack < 0:
                        print("COUNTEREXAMPLE", n, restart, step, cur.cap_slack, cur.q_slack)
                        print(A.tolist())
                        print(cur.data[:-1])
                        return 1
            print("hill-final", n, "cap", best[0], "q", best[1])
        return 0
    global_q = (10**18, None)
    global_cap = (10**18, None)
    for n in range(args.n_min, args.n_max + 1):
        nq = (10**18, None)
        nc = (10**18, None)
        for trial in range(args.trials):
            A = random_signing(n, rng)
            au = audit_matrix(A, cache)
            if au.q_slack < nq[0]:
                nq = (au.q_slack, (trial, au.data))
            if au.cap_slack < nc[0]:
                nc = (au.cap_slack, (trial, au.data))
            if au.q_slack < 0 or au.cap_slack < 0:
                print("COUNTEREXAMPLE", n, trial, au.q_slack, au.cap_slack)
                print(A.tolist())
                print(au.data[:-1])
                return 1
        global_q = min(global_q, (nq[0], (n, nq[1])), key=lambda z: z[0])
        global_cap = min(global_cap, (nc[0], (n, nc[1])), key=lambda z: z[0])
        print("n", n, "min_q_slack", nq[0], "min_cap_slack", nc[0])
    print("global_q", global_q[0], "at", global_q[1][0:1])
    print("global_cap", global_cap[0], "at", global_cap[1][0:1])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
