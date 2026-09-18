#!/usr/bin/env python3
"""Exploratory circle-valued low-degree least squares, NOT a certificate."""
import argparse
import json
import math
import numpy as np
from scipy.optimize import least_squares


def weight(x):
    return bin(int(x)).count("1")


def run(n, m, starts, seed, balanced=False):
    rng = np.random.RandomState(seed)
    masks = [s for s in range(1 << n) if weight(s) <= m and (s or not balanced)]
    W = np.array([[(-1) ** weight(x & s) for s in masks]
                  for x in range(1 << n)], dtype=float)
    d = len(masks)
    rows = []
    best = None
    best_feasible = None
    for trial in range(starts):
        ar = rng.normal(size=d)
        ai = rng.normal(size=d)
        scale = math.sqrt(float(ar @ ar + ai @ ai))
        z0 = np.r_[ar, ai] / scale

        def fun(z):
            fr, fi = W @ z[:d], W @ z[d:]
            return fr * fr + fi * fi - 1

        def jac(z):
            fr, fi = W @ z[:d], W @ z[d:]
            return np.concatenate((2 * fr[:, None] * W,
                                   2 * fi[:, None] * W), axis=1)

        opt = least_squares(fun, z0, jac=jac, gtol=1e-11, ftol=1e-12,
                            xtol=1e-12, max_nfev=1200)
        a = opt.x[:d] + 1j * opt.x[d:]
        q = 2 * m / (m + 1.0)
        supremum = float(np.max(np.abs(W @ a)))
        ratio = float(np.sum(np.abs(a) ** q) ** (1 / q) / supremum)
        residual = float(np.max(np.abs(fun(opt.x))))
        support = int(np.count_nonzero(np.abs(a) > 1e-6))
        row = {"trial": trial, "residual": residual, "support_1e-6": support,
               "ratio": ratio, "nfev": opt.nfev}
        rows.append(row)
        print(json.dumps({"n": n, "m": m, **row}), flush=True)
        if best is None or ratio > best[0]:
            best = (ratio, a, residual)
        if residual < 1e-8 and (best_feasible is None or ratio > best_feasible[0]):
            best_feasible = (ratio, a, residual)
    print(json.dumps({"summary": {"n": n, "m": m, "starts": starts,
          "best_ratio": best[0], "best_residual": best[2],
          "coefficients": [[s, float(a.real), float(a.imag)]
                           for s, a in zip(masks, best[1])]}}), flush=True)
    print(json.dumps({"feasibility_summary": {
          "criterion": "max modulus-squared residual < 1e-8; numerical only",
          "count": sum(row["residual"] < 1e-8 for row in rows),
          "best_ratio": None if best_feasible is None else best_feasible[0],
          "best_residual": None if best_feasible is None else best_feasible[2],
          "largest_support": max([row["support_1e-6"] for row in rows
                                   if row["residual"] < 1e-8] + [0])}}), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=6)
    parser.add_argument("--m", type=int, default=2)
    parser.add_argument("--starts", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20260918)
    parser.add_argument("--balanced", action="store_true")
    args = parser.parse_args()
    run(args.n, args.m, args.starts, args.seed, args.balanced)
