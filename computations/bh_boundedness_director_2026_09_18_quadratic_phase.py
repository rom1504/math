#!/usr/bin/env python3
"""Diagnostic only: all-time quadratic phase entropy and rotation speed.

The theorem is proved independently in director artifact Section 7.
This replay deliberately includes gauge shifts and zero cosine factors.
"""
import argparse
import json
from pathlib import Path

import numpy as np


def fwht(values):
    out = np.asarray(values, dtype=complex).copy()
    step = 1
    while step < len(out):
        blocks = out.reshape(-1, 2 * step)
        left, right = blocks[:, :step].copy(), blocks[:, step:].copy()
        blocks[:, :step], blocks[:, step:] = left + right, left - right
        step *= 2
    return out / len(out)


def entropy(p):
    p = np.asarray(p)
    nz = p > 1e-28
    return float(-np.sum(p[nz] * np.log(p[nz])))


def check(theta, linear=None):
    n = len(theta)
    if linear is None:
        linear = np.zeros(n)
    vertices = np.arange(1 << n)
    spins = 1 - 2 * ((vertices[:, None] >> np.arange(n)) & 1)
    energy = np.einsum("bi,ij,bj->b", spins, theta, spins) / 2 + spins @ linear
    g = np.exp(1j * energy)
    prob = np.abs(fwht(g)) ** 2
    H = entropy(prob)
    J = np.zeros((1 << n, 1 << n), complex)
    for i in range(n):
        J[vertices, vertices ^ (1 << i)] = 0.5
    A = g.conj()[:, None] * J * g[None, :] - J
    v = float(np.max(np.abs(np.linalg.eigvalsh(A))))
    residual = (theta + np.pi / 4) % (np.pi / 2) - np.pi / 4
    np.fill_diagonal(residual, 0)
    residual_linear = (linear + np.pi/4) % (np.pi/2) - np.pi/4
    c = np.cos(2 * residual_linear) * np.prod(np.cos(2 * residual), axis=1)
    q = np.maximum(0, (1 - c) / 2)
    hollow = bool(np.all(linear == 0))
    small = q <= (0.125 if hollow else 1/32) + 1e-13
    small_witness = (3/64 if hollow else 1/16) * np.sum(np.sqrt(q[small]))
    marg = np.array([prob[((vertices >> i) & 1) == 1].sum()
                     for i in range(n)])
    marginal_entropy = sum(entropy([p, 1-p]) for p in marg)
    assert H <= marginal_entropy + 1e-9
    assert v + 1e-9 >= q.sum()
    assert v + 1e-9 >= small_witness
    assert H <= (152 / 3 if hollow else 64) * v + 1e-9
    # Check residual reduction preserves marginal entropy, not velocity.
    assert abs(sum(entropy([p, 1-p]) for p in q)
               - marginal_entropy) < 1e-8
    return {"n": n, "theta": theta.tolist(), "linear": linear.tolist(), "entropy": H,
            "velocity": v, "entropy_over_velocity": H/v if v > 1e-12 else None,
            "marginal_entropy": marginal_entropy,
            "sum_residual_marginals": float(q.sum()),
            "small_row_witness": float(small_witness)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rng = np.random.default_rng(9182026)
    rows = []
    for n in range(2, 7):
        for scale in (0, 1e-4, 0.03, 0.2, 1, 4, 30):
            for repeat in range(3):
                theta = rng.normal(size=(n, n)) * scale
                theta = np.triu(theta, 1)
                theta += theta.T
                rows.append(check(theta))
                rows.append(check(theta, rng.normal(size=n)*scale))
        for angle in (np.pi/4, np.pi/2, 3*np.pi/4):
            theta = np.full((n, n), angle)
            np.fill_diagonal(theta, 0)
            rows.append(check(theta))
    result = {"status": "NUMERICAL_DIAGNOSTIC_ONLY", "seed": 9182026,
              "checks": len(rows), "rows": rows}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"checks": len(rows), "status": "PASS_DIAGNOSTICS",
                      "maximum_observed_ratio": max(
                          r["entropy_over_velocity"] or 0 for r in rows)}))


if __name__ == "__main__":
    main()
