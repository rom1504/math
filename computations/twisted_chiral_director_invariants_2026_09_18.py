#!/usr/bin/env python3
"""Exact twist-invariant bilinear and diagonal-optimized necessary bounds.

Every integer entry in the output is exhaustive over the stated finite cube.
These are necessary bounds for a twisted double, not sufficient cap estimates.
"""
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def invariant(a):
    a = np.asarray(a, dtype=np.int64)
    n = len(a)
    spins = 1 - 2 * ((np.arange(1 << n)[:, None] >> np.arange(n)) & 1)
    local = spins * (spins @ a)
    l1 = np.abs(local).sum(axis=1)
    slopes = np.sign(local)
    zero = (local == 0).sum(axis=1)
    # |integer h + sign d| = |h| + 1_{h=0} + sign(h)d.
    values = l1[:, None] + zero[:, None] + slopes @ spins.T
    caps = values.max(axis=0)
    best = int(caps.argmin())
    d = spins[best]
    direct = np.abs(spins @ (a + np.diag(d))).sum(axis=1)
    assert np.array_equal(direct, values[:, best])
    q = int(np.abs(local.sum(axis=1) // 2).max())
    hist = {str(int(k)): int(v) for k, v in zip(*np.unique(caps, return_counts=True))}
    return {"order": n, "child_cap": q, "bilinear_norm": int(l1.max()),
            "diagonal_optimized_bilinear_norm": int(caps.min()),
            "diagonal_cap_histogram": hist, "best_diagonal": d.tolist(),
            "active_bilinear_state_count": int((l1 == l1.max()).sum()),
            "matrix": a.tolist()}


def main():
    records = []
    for n in range(3, 11):
        source = ROOT / f"computations/results/m{n}_minimizer_orbits.json"
        if source.exists():
            data = json.loads(source.read_text())
            seeds = [(c["class"], c["representative_matrix"]) for c in data["classes"]]
        else:
            source = ROOT / f"computations/results/exact_m{n}.json"
            seeds = [(0, json.loads(source.read_text())["matrix"])]
        for label, a in seeds:
            r = invariant(a)
            r.update(source=str(source.relative_to(ROOT)), seed_class=label,
                     evidence="exact exhaustive finite necessary bound; no global/family-optimum inference")
            records.append(r)
    print(json.dumps({"records": records}, indent=2))


if __name__ == "__main__":
    main()
