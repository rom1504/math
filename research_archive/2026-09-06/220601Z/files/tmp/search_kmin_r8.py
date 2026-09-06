#!/usr/bin/env python3
"""Switching-gauge and random search for large optimized K_min."""

import argparse
from multiprocessing import Pool
from random import Random

from kmin_carleson_r8 import EndpointSystem, Model


def gauge_matrix(n, code):
    A = [[0] * n for _ in range(n)]
    for j in range(1, n):
        A[0][j] = A[j][0] = 1
    bit = 0
    for i in range(1, n):
        for j in range(i + 1, n):
            v = 1 if code >> bit & 1 else -1
            A[i][j] = A[j][i] = v
            bit += 1
    return tuple(tuple(row) for row in A)


def solve_one(arg):
    n, code, seconds = arg
    A = gauge_matrix(n, code)
    S = EndpointSystem(A)
    M = Model(S)
    res = M.solve(time_limit=seconds, mip_gap=1e-9)
    return code, res.status, res.fun, getattr(res, "mip_gap", None), A


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=6)
    ap.add_argument("--trials", type=int, default=0, help="0 means exhaustive switching gauge")
    ap.add_argument("--seed", type=int, default=81623)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--seconds", type=float, default=30)
    args = ap.parse_args()
    total = 1 << ((args.n - 1) * (args.n - 2) // 2)
    if args.trials:
        rng = Random(args.seed)
        codes = [rng.randrange(total) for _ in range(args.trials)]
    else:
        codes = range(total)
    best = (-1, None)
    solved = 0
    with Pool(args.jobs) as pool:
        for code, status, value, gap, A in pool.imap_unordered(
            solve_one, ((args.n, c, args.seconds) for c in codes), chunksize=1
        ):
            if status == 0:
                solved += 1
                if value > best[0] + 1e-9:
                    best = (value, code, A)
                    print("BEST", args.n, solved, code, value, flush=True)
                    for row in A:
                        print(row, flush=True)
            else:
                print("UNSOLVED", code, status, value, gap, flush=True)
    print("DONE", args.n, "solved", solved, "best", best[:2])


if __name__ == "__main__":
    main()
