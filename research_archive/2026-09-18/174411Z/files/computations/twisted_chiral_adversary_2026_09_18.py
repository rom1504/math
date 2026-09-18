#!/usr/bin/env python3
"""Independent exact checks for the 2026-09-18 twisted-chiral audit.

Uses H=sum_{i<j}; no heuristic cap is reported as an exact cap.
"""

import argparse
import json
from pathlib import Path

import numpy as np


def spins(n):
    bits = (np.arange(1 << n, dtype=np.uint64)[:, None]
            >> np.arange(n, dtype=np.uint64)[None, :]) & 1
    return (2 * bits.astype(np.int64) - 1)


def energies(a, x):
    return np.einsum("bi,ij,bj->b", x, a, x, optimize=True) // 2


def direct_cap(a):
    x = spins(len(a))
    h = energies(a, x)
    i = int(np.argmax(np.abs(h)))
    return int(abs(h[i])), x[i].tolist(), int(h[i])


def partition_cap(a, b, d):
    n = len(a)
    z = spins(n)
    best = -1
    witness = None
    for mask in range(1 << n):
        inside = ((mask >> np.arange(n)) & 1).astype(bool)
        outside = ~inside
        cross = np.einsum("bi,ij,bj->b", z[:, inside],
                          a[np.ix_(inside, outside)], z[:, outside],
                          optimize=True)
        inner = (2 * energies(b[np.ix_(inside, inside)], z[:, inside])
                 - 2 * energies(b[np.ix_(outside, outside)], z[:, outside])
                 + int(d[inside].sum() - d[outside].sum()))
        vals = 2 * abs(cross) + abs(inner)
        k = int(np.argmax(vals))
        if vals[k] > best:
            best = int(vals[k])
            witness = {"I": np.flatnonzero(inside).tolist(),
                       "z": z[k].tolist(), "cut": int(cross[k]),
                       "remainder": int(inner[k])}
    return best, witness


def random_signing(n, rng):
    a = rng.choice([-1, 1], size=(n, n))
    a = np.triu(a, 1)
    return a + a.T


def run_checks(seed=202609181711):
    rng = np.random.default_rng(seed)
    checks = []
    for n in range(1, 7):
        for case in range(8):
            a = random_signing(n, rng)
            b = random_signing(n, rng)
            d = rng.choice([-1, 1], size=n)
            c = b + np.diag(d)
            double = np.block([[a, c], [c, -a]])
            exact, _, _ = direct_cap(double)
            formula, witness = partition_cap(a, b, d)
            assert exact == formula
            checks.append({"n": n, "case": case, "cap": exact,
                           "A": a.tolist(), "B": b.tolist(),
                           "d": d.tolist(), "formula_witness": witness})

    # Integer-scaled exterior-escape example, denominator 10.
    w = np.zeros((4, 4), dtype=np.int64)
    w[0, 1:] = w[1:, 0] = 1
    v = np.zeros_like(w)
    v[1, 2] = v[2, 1] = 10
    v[1, 3] = v[3, 1] = -10
    x = spins(4)
    old = energies(w, x)
    change = energies(v, x)
    assert int(abs(old).max()) == 3
    assert np.all(change[abs(old) == 3] == 0)
    assert int(abs(old + change).max()) == 21

    # Exact edge-coordinate mixed differences.
    shell_checks = []
    for n in range(2, 15):
        rows, cols = np.triu_indices(n, 1)
        anchor = rng.choice([-1, 1], size=n)
        feature = lambda z: z[rows] * z[cols]
        base = feature(anchor)
        for edge, (i, j) in enumerate(zip(rows, cols)):
            one, two, both = anchor.copy(), anchor.copy(), anchor.copy()
            one[i] *= -1
            two[j] *= -1
            both[[i, j]] *= -1
            diff = base - feature(one) - feature(two) + feature(both)
            expected = np.zeros(len(rows), dtype=np.int64)
            expected[edge] = 4 * anchor[i] * anchor[j]
            assert np.array_equal(diff, expected)
        shell_checks.append({"n": n, "edge_rank": len(rows)})
    rounding_checks = []
    for m in range(3, 6):
        n = 3*m
        group = np.arange(n)//m
        a = np.where(group[:, None] == group[None, :], -1, 1)
        np.fill_diagonal(a, 0)
        denominator = 3*m-1
        wnum = (m+1)*(np.ones((n, n), dtype=np.int64)-np.eye(n, dtype=np.int64))
        vnum = denominator*a-wnum
        assert np.all(vnum.sum(axis=1) == 0)
        xx = spins(n)
        hw = energies(wnum, xx)
        hv = energies(vnum, xx)
        ha = energies(a, xx)
        flip = (1-xx)//2
        assert np.array_equal(hv, 2*np.einsum("bi,ij,bj->b", flip, vnum, flip, optimize=True))
        assert int(abs(ha).max()) == (5*m*m-3*m)//2
        rounding_checks.append({"m": m, "n": n, "rho_numerator": m+1,
                                "rho_denominator": denominator,
                                "old_cap": (3*m*m+3*m)//2,
                                "rounded_cap": int(abs(ha).max()),
                                "row_sum_zero_verified": True,
                                "flip_identity_verified_states": len(xx)})
    return {"seed": seed, "normalization": "sum_i<j", "identity_checks": checks,
            "shell_mixed_difference_checks": shell_checks,
            "approximate_shell_full_rounding_checks": rounding_checks,
            "escape_example": {"denominator": 10, "W": w.tolist(),
                               "V": v.tolist(), "old_cap_numerator": 3,
                               "new_cap_numerator": 21}}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run_checks()
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"identity_checks": len(result["identity_checks"]),
                      "shell_orders": len(result["shell_mixed_difference_checks"]),
                      "escape_old_cap": "3/10", "escape_new_cap": "21/10"}))


if __name__ == "__main__":
    main()
