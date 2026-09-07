#!/usr/bin/env python3
"""Exploratory finite soft-center capacities for A6/A8/A9.

Numerical optimization is evidence only.  Geometry and reported candidate
overlaps are enumerated exactly before exponentiation.
"""

from __future__ import annotations

import itertools
import math
import sys

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, "/home/math/quadra/tmp")
from envelope_block_cover_r27 import A6, A8, A9, projective_spins


def pspins(n):
    return np.array(list(projective_spins(n)), dtype=np.int8)


def pdist(x, y):
    d = int(np.sum(x != y))
    return min(d, len(x) - d)


def geometry(A, m, tolerance=0.0):
    n = len(A)
    zs = pspins(n)
    row = np.einsum("bi,ij,bj->b", zs, A @ A, zs)
    centers = zs[row <= 2 * n * (n - 1)]
    selectors = list(itertools.combinations(range(n), m))
    q = int(np.max(np.abs(np.einsum("bi,ij,bj->b", zs, A, zs))))
    p = m / n
    p2 = m * (m - 1) / (n * (n - 1))
    B = (p ** 1.5 - p2) * q
    distances = np.empty((len(centers), len(selectors)), dtype=np.int8)
    fav_counts = []
    qS_list = []
    for j, S in enumerate(selectors):
        S = list(S)
        ys = pspins(m)
        AS = A[np.ix_(S, S)]
        en = np.abs(np.einsum("bi,ij,bj->b", ys, AS, ys))
        qS = int(en.max())
        fav = ys[(qS - en) <= B + tolerance + 1e-12]
        qS_list.append(qS)
        fav_counts.append(len(fav))
        for i, z in enumerate(centers):
            distances[i, j] = min(pdist(z[S], y) for y in fav)
    return {
        "centers": centers,
        "selectors": selectors,
        "dist": distances,
        "row": row[row <= 2 * n * (n - 1)],
        "q": q,
        "B": B,
        "fav_counts": fav_counts,
        "qS": qS_list,
    }


def capacity(K, s, starts=6):
    nz, ns = K.shape
    nu = 1.0 / nz
    def fun(w):
        u = K @ w
        moment = nu * np.sum(u ** s)
        return moment ** (1.0 / s)
    def jac(w):
        u = K @ w
        moment = nu * np.sum(u ** s)
        # derivative of (E u^s)^(1/s)
        return (moment ** (1.0 / s - 1.0)) * nu * (K.T @ (u ** (s - 1)))
    candidates = [np.full(ns, 1 / ns)]
    rng = np.random.default_rng(3402)
    for _ in range(starts - 1):
        candidates.append(rng.dirichlet(np.ones(ns)))
    best = None
    for x0 in candidates:
        res = minimize(fun, x0, jac=jac, method="SLSQP",
                       constraints={"type": "eq", "fun": lambda w: np.sum(w)-1,
                                    "jac": lambda w: np.ones_like(w)},
                       bounds=[(0, 1)] * ns,
                       options={"ftol": 2e-12, "maxiter": 3000})
        if best is None or res.fun < best.fun:
            best = res
    w = best.x
    u = K @ w
    V = fun(w)
    h = u ** (s - 1) / V ** (s - 1)
    overlaps = (h[:, None] * K).mean(axis=0)
    return best, w, h, overlaps


def slack_candidate(A, geo, lam, theta, s):
    z = geo["centers"]
    en = np.abs(np.einsum("bi,ij,bj->b", z, A, z))
    delta = geo["q"] - en
    raw = np.exp(-theta * delta)
    qexp = s / (s - 1)
    h = raw / (np.mean(raw ** qexp) ** (1 / qexp))
    K = np.exp(-lam * geo["dist"])
    return float(np.min(np.mean(h[:, None] * K, axis=0))), h


def run(name, A, m):
    geo = geometry(A, m)
    d = geo["dist"]
    print("GEOM", name, "n,m", len(A), m, "q", geo["q"], "B", geo["B"],
          "centers", len(geo["centers"]), "selectors", len(geo["selectors"]),
          "dmax", int(d.max()), "fav", (min(geo["fav_counts"]), max(geo["fav_counts"])),
          "qS", sorted(set(geo["qS"])))
    for s in (2, 3, 5, 8):
        for lam in (0.25, 0.5, 1.0, 2.0, 4.0):
            K = np.exp(-lam * d)
            res, w, h, ov = capacity(K, s)
            slack = max(slack_candidate(A, geo, lam, theta, s)[0]
                        for theta in (0, .05, .1, .2, .5, 1, 2, 4))
            print("CAP", name, m, "s", s, "lam", lam,
                  "V", f"{res.fun:.12g}", "minov", f"{ov.min():.12g}",
                  "maxov", f"{ov.max():.12g}", "support", int(np.sum(w > 1e-7)),
                  "slackbest", f"{slack:.12g}", "ok", res.success,
                  "nit", res.nit)


if __name__ == "__main__":
    mats = {"A6": A6, "A8": A8, "A9": A9}
    if len(sys.argv) != 3:
        raise SystemExit("usage: script A6|A8|A9 m")
    run(sys.argv[1], mats[sys.argv[1]], int(sys.argv[2]))
