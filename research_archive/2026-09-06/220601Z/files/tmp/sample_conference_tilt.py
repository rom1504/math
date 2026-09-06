#!/usr/bin/env python3
"""Seeded Monte Carlo for the conference-child tilted bridge landscape."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def spins(n):
    grid = np.arange(1 << n, dtype=np.uint64)[:, None]
    bits = (grid >> np.arange(n, dtype=np.uint64)) & 1
    return (2 * bits.astype(np.int8) - 1)


def energy(a, x):
    return np.einsum("bi,ij,bj->b", x, a, x, dtype=np.int64) // 2


def logmean_cosh(v):
    v = np.abs(np.asarray(v, dtype=float))
    p = float(np.max(v))
    return p + math.log(float(np.mean((np.exp(v - p) + np.exp(-v - p)) / 2)))


def logmeanexp(v):
    p = float(np.max(v))
    return p + math.log(float(np.mean(np.exp(v - p))))


def main():
    r = 6
    samples = 20_000
    rng = np.random.default_rng(20260815)
    payload = json.loads(
        (ROOT / "computations/results/conference_double_p5.json").read_text()
    )
    a = np.asarray(payload["conference_matrix"], dtype=np.int8)
    x = spins(r)
    e = energy(a, x)
    rows = []
    for beta in (0.2, 0.5, 1.0):
        t = beta / math.sqrt(2 * r)
        s = beta / math.sqrt(r)
        target = 2 * logmean_cosh(s * e)
        vals = np.empty(2 * samples)
        for k in range(samples):
            b = rng.choice(np.asarray([-1, 1], dtype=np.int8), size=(r, r))
            cross = x.astype(np.int64) @ b.astype(np.int64) @ x.astype(np.int64).T
            vals[2 * k] = logmean_cosh(t * (e[:, None] - e[None, :] + cross).ravel())
            vals[2 * k + 1] = logmean_cosh(t * (e[:, None] + e[None, :] + cross).ravel())
        rows.append({
            "beta": beta,
            "samples": samples,
            "target": target,
            "mean_defect": float(np.mean(vals) - target),
            "minimum_sample_defect": float(np.min(vals) - target),
            "quantiles_defect": {
                str(q): float(np.quantile(vals, q) - target)
                for q in (0.001, 0.01, 0.1, 0.5, 0.9)
            },
            "tilted_defects": {
                str(lam): float(-logmeanexp(-lam * vals) / lam - target)
                for lam in (0.25, 1.0, 4.0, 16.0, 64.0)
            },
        })
    print(json.dumps({"classification": "seeded Monte Carlo", "rows": rows}, indent=2))


if __name__ == "__main__":
    main()
