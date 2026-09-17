#!/usr/bin/env python3
"""Exhaustive finite-order audit of entropy-tilted bridge soft minima.

This is a scratch diagnostic.  It enumerates switching-normalized child
signings, all bridge bits, and both relative orientations.  The output is
floating-point, while the finite enumerations are exhaustive.
"""

from __future__ import annotations

import itertools
import json
import math

import numpy as np


def spins(n: int) -> np.ndarray:
    return np.asarray(list(itertools.product((-1, 1), repeat=n)), dtype=np.int8)


def signings(n: int) -> list[np.ndarray]:
    if n <= 1:
        return [np.zeros((n, n), dtype=np.int8)]
    edges = [(i, j) for i in range(1, n) for j in range(i + 1, n)]
    out = []
    for mask in range(1 << len(edges)):
        a = np.zeros((n, n), dtype=np.int8)
        a[0, 1:] = a[1:, 0] = 1
        for bit, (i, j) in enumerate(edges):
            a[i, j] = a[j, i] = -1 if (mask >> bit) & 1 else 1
        out.append(a)
    return out


def energy(a: np.ndarray, x: np.ndarray) -> np.ndarray:
    return np.einsum("bi,ij,bj->b", x, a, x, dtype=np.int64) // 2


def logmean_cosh(v: np.ndarray) -> float:
    v = np.abs(np.asarray(v, dtype=float))
    p = float(np.max(v))
    return p + math.log(float(np.mean((np.exp(v - p) + np.exp(-v - p)) / 2)))


def pressure(e: np.ndarray, t: float) -> float:
    return logmean_cosh(t * e)


def logmeanexp(v: np.ndarray) -> float:
    p = float(np.max(v))
    return p + math.log(float(np.mean(np.exp(v - p))))


def landscape(ae, de, x, y, t):
    m, n = x.shape[1], y.shape[1]
    ia, id_ = ae[:, None], de[None, :]
    vals = []
    for mask in range(1 << (m * n)):
        flat = np.ones(m * n, dtype=np.int8)
        for bit in range(m * n):
            if (mask >> bit) & 1:
                flat[bit] = -1
        b = flat.reshape(m, n)
        cross = x.astype(np.int64) @ b.astype(np.int64) @ y.astype(np.int64).T
        for eps in (-1, 1):
            vals.append(pressure((ia + eps * id_ + cross).ravel(), t))
    return np.asarray(vals)


def main():
    max_n = 7
    betas = (0.25, 0.5, 1.0, 2.0)
    lambdas = (0.25, 1.0, 4.0, 16.0)
    xs = {n: spins(n) for n in range(1, max_n + 1)}
    mats = {n: signings(n) for n in range(1, max_n + 1)}
    es = {n: [energy(a, xs[n]) for a in mats[n]] for n in mats}
    records = []
    for beta in betas:
        high = {}
        for n in range(1, max_n + 1):
            vals = np.asarray([pressure(e, beta / math.sqrt(n)) for e in es[n]])
            high[n] = float(np.min(vals))
        for total in range(2, max_n + 1):
            m = total // 2
            n = total - m
            t = beta / math.sqrt(total)
            mv = np.asarray([pressure(e, t) for e in es[m]])
            nv = np.asarray([pressure(e, t) for e in es[n]])
            mi = np.flatnonzero(mv <= np.min(mv) + 1e-11)
            ni = np.flatnonzero(nv <= np.min(nv) + 1e-11)
            best = None
            for i in mi:
                for j in ni:
                    ls = landscape(es[m][i], es[n][j], xs[m], xs[n], t)
                    row = {
                        "mean_log": float(np.mean(ls)),
                        "minimum_log": float(np.min(ls)),
                        "tilted": {
                            str(lam): float(-logmeanexp(-lam * ls) / lam)
                            for lam in lambdas
                        },
                    }
                    score = tuple(row["tilted"][str(lam)] for lam in lambdas)
                    if best is None or score < best[0]:
                        best = (score, row)
            row = best[1]
            target = high[m] + high[n]
            records.append({
                "beta": beta,
                "N": total,
                "split": [m, n],
                "target_same_beta": target,
                "mean_defect": row["mean_log"] - target,
                "best_bridge_defect": row["minimum_log"] - target,
                "tilted_defects": {
                    key: value - target for key, value in row["tilted"].items()
                },
            })
    print(json.dumps({
        "classification": "exhaustive enumeration; floating-point evaluation",
        "records": records,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
