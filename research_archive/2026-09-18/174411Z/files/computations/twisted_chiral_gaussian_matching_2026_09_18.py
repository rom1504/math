#!/usr/bin/env python3
"""Exact coefficient/covariance audit for signed-complex-structure twists.

No Gaussian random sampling is needed: covariance is a coefficient dot product.
The output tests the natural comparison to two independent copies of the child
Gaussian process, and retains explicit witnesses when covariance order fails.
"""

import itertools
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def spins(n):
    return np.asarray(list(itertools.product((-1, 1), repeat=n)), dtype=np.int64)


def gaussian_coefficients(x, y, g):
    u, v = g @ x, g @ y
    return np.asarray([
        x[i] * x[j] - y[i] * y[j] + u[i] * v[j] + u[j] * v[i]
        for i in range(len(x)) for j in range(i)
    ], dtype=np.int64)


def parent_covariance(x, y, xp, yp, g):
    a, b, c, d = x @ xp, x @ yp, y @ xp, y @ yp
    ag, bg, cg, dg = x @ g @ xp, x @ g @ yp, y @ g @ xp, y @ g @ yp
    r = (x * y) @ (xp * yp)
    return int(((a + d) ** 2 - (b - c) ** 2) // 2 + (ag - dg) * (bg + cg) - 2 * r)


def reference_covariance(x, y, xp, yp):
    return int((x @ xp) ** 2 + (y @ yp) ** 2 - 2 * len(x))


def iid_core_covariance(x, y, xp, yp):
    n = len(x)
    return int(((x @ xp) + (y @ yp)) ** 2 // 2 - n - (x*y) @ (xp*yp))


def main():
    n = 4
    g = np.asarray([[0, -1, 0, 0], [1, 0, 0, 0],
                    [0, 0, 0, -1], [0, 0, 1, 0]], dtype=np.int64)
    assert np.array_equal(g.T, -g)
    assert np.array_equal(g @ g, -np.eye(n, dtype=np.int64))
    cube = spins(n)
    states = [(x, y) for x in cube for y in cube]
    coeffs = np.asarray([gaussian_coefficients(x, y, g) for x, y in states])
    parent_cov = coeffs @ coeffs.T
    reference_cov = np.asarray([
        [reference_covariance(x, y, xp, yp) for xp, yp in states]
        for x, y in states
    ], dtype=np.int64)
    for i, (x, y) in enumerate(states):
        assert parent_cov[i, i] == 2 * n * (n - 1)
        for j, (xp, yp) in enumerate(states):
            assert parent_covariance(x, y, xp, yp, g) == parent_cov[i, j]
    differences = parent_cov - reference_cov
    witnesses = []
    for kind, flat in (("negative", differences.argmin()), ("positive", differences.argmax())):
        i, j = np.unravel_index(flat, differences.shape)
        x, y = states[i]
        xp, yp = states[j]
        replications = []
        for copies in (1, 2, 4, 8, 16):
            gg = np.kron(np.eye(copies, dtype=np.int64), g)
            xx, yy, xxp, yyp = [np.tile(z, copies) for z in (x, y, xp, yp)]
            p = parent_covariance(xx, yy, xxp, yyp, gg)
            r = reference_covariance(xx, yy, xxp, yyp)
            replications.append({"n": len(xx), "parent_covariance": p,
                                 "reference_covariance": r, "difference": p-r})
        witnesses.append({"kind": kind, "x": x.tolist(), "y": y.tolist(),
                          "xp": xp.tolist(), "yp": yp.tolist(),
                          "parent_coefficients": coeffs[i].tolist(),
                          "parent_prime_coefficients": coeffs[j].tolist(),
                          "replications": replications})
    result = {
        "description": "Equal-variance matching-twist Gaussian process versus sqrt(2) times the sum of two independent child processes",
        "n": n, "g": g.tolist(), "states": len(states),
        "covariance_entries_checked": int(parent_cov.size),
        "common_variance": 2*n*(n-1),
        "minimum_covariance_difference": int(differences.min()),
        "maximum_covariance_difference": int(differences.max()),
        "witnesses": witnesses,
    }
    iid_covariance = np.asarray([
        [iid_core_covariance(x, y, xp, yp) for xp, yp in states]
        for x, y in states
    ], dtype=np.int64)
    iid_differences = parent_cov - iid_covariance
    balanced = [i for i, (x, y) in enumerate(states) if x @ y == 0]
    iid_witnesses = []
    for kind, direction in (("negative", -1), ("positive", 1)):
        candidates = []
        for i in balanced:
            for j in balanced:
                x, y = states[i]
                xp, yp = states[j]
                # Exclude the trivial global/chiral symmetry orbit, so the
                # increment failure is not only a redundant antipodal test.
                chiral = ((x, y), (y, -x), (-x, -y), (-y, x))
                if any(np.array_equal(xp, xx) and np.array_equal(yp, yy) for xx, yy in chiral):
                    continue
                candidates.append((direction * iid_differences[i, j], i, j))
        _, i, j = max(candidates)
        x, y = states[i]
        xp, yp = states[j]
        replications = []
        for copies in (1, 2, 4, 8, 16):
            gg = np.kron(np.eye(copies, dtype=np.int64), g)
            xx, yy, xxp, yyp = [np.tile(z, copies) for z in (x, y, xp, yp)]
            p = parent_covariance(xx, yy, xxp, yyp, gg)
            r = iid_core_covariance(xx, yy, xxp, yyp)
            replications.append({"n": len(xx), "parent_covariance": p,
                                 "reference_covariance": r, "difference": p-r})
        iid_witnesses.append({"kind": kind, "x": x.tolist(), "y": y.tolist(),
                              "xp": xp.tolist(), "yp": yp.tolist(),
                              "replications": replications})
    result["iid_core_comparison"] = {
        "minimum_covariance_difference": int(iid_differences.min()),
        "maximum_covariance_difference": int(iid_differences.max()),
        "balanced_nonchiral_witnesses": iid_witnesses,
    }
    output = ROOT / "computations/results/twisted_chiral_gaussian_matching_2026_09_18.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "witnesses"}))
    for witness in witnesses:
        print(witness["kind"], witness["x"], witness["y"], witness["xp"], witness["yp"], witness["replications"])
    for witness in iid_witnesses:
        print("iid_core_balanced", witness)


if __name__ == "__main__":
    main()
