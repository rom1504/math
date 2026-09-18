#!/usr/bin/env python3
"""Numerical rotation-speed diagnostics for the exact ten-bit tangent seed.

Numerical eigenvalues are diagnostics, not proof certificates. Exact seed
algebra is independently checked by the companion tangent replay.
"""
import argparse
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import eigsh


def load_exact_helpers():
    path = Path(__file__).with_name("bh_boundedness_barrier_2026_09_18_tangent.py")
    spec = importlib.util.spec_from_file_location("exact_tangent", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def adjacency(n):
    vertices = np.arange(1 << n)
    rows = np.tile(vertices, n)
    cols = np.concatenate([vertices ^ (1 << j) for j in range(n)])
    return coo_matrix((np.full(len(rows), 0.5), (rows, cols)),
                      shape=(1 << n, 1 << n)).tocsr()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    t = load_exact_helpers()
    roots = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
    masks = [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 16, 17, 18, 20, 24]
    powers = [0, 2, 0, 5, 1, 4, 2, 4, 3, 1, 4, 1, 2, 2, 3, 2]
    f4 = dict(zip(masks, [roots[j] for j in powers]))
    a = {0: (1, -1), 4: (0, 1), 8: (1, 0)}
    b16 = t.convolution(t.convolution(f4, f4),
                        {i: t.conj(z) for i, z in a.items()})
    F16 = t.tensor(f4, f4)
    h16 = t.tensor(a, b16)
    for mask, z in t.tensor(b16, a).items():
        h16[mask] = t.add(h16.get(mask, (0, 0)), t.neg(z))
    zeta = np.exp(1j*np.pi/3)

    def values(poly, denominator, n):
        pairs = [t.evaluate(poly, x) for x in range(1 << n)]
        return np.array([u+v*zeta for u, v in pairs])/denominator

    fv = values(f4, 4, 5)
    Fv, hv = values(F16, 16, 10), values(h16, 16, 10)
    J5, J10 = adjacency(5), adjacency(10)
    A5 = (diags(fv.conj()) @ J5 @ diags(fv)-J5).toarray()
    e5, u5 = np.linalg.eigh(A5)
    peak5 = float(e5[-1])
    peak_basis5 = u5[:, e5 > peak5-1e-9]
    peak_basis10 = np.kron(peak_basis5, peak_basis5)
    phase_tangent = np.real(-1j*hv/Fv)
    B = diags(Fv.conj()) @ J10 @ diags(Fv)
    T = diags(phase_tangent)
    derivative = 1j*(B @ T-T @ B)
    compression = peak_basis10.conj().T @ (derivative @ peak_basis10)
    rows = []
    for epsilon in (0.0, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 0.03, 0.1):
        P = Fv+epsilon*hv
        G = P/np.abs(P)
        A = diags(G.conj()) @ J10 @ diags(G)-J10
        vals, vectors = eigsh(A, k=4, which="LA", tol=2e-13,
                              maxiter=5000, v0=np.ones(1024))
        order = np.argsort(vals)
        vals, vectors = vals[order], vectors[:, order]
        residual = np.linalg.norm(A @ vectors[:, -1]-vals[-1]*vectors[:, -1])
        row = {"epsilon": epsilon, "largest_eigenvalues": vals.tolist(),
               "eigen_residual": float(residual),
               "velocity": float(vals[-1]),
               "difference_from_base": float(vals[-1]-2*peak5),
               "quadratic_quotient": (float((vals[-1]-2*peak5)/epsilon**2)
                                       if epsilon else None),
               "effective_degree_upper_bound":
                    4+8*np.log1p(48*epsilon**2)}
        rows.append(row)
        print(json.dumps(row), flush=True)
    result = {"status": "NUMERICAL_DIAGNOSTIC_ONLY",
              "base_f_velocity": peak5,
              "base_f_spectrum": e5.tolist(),
              "base_top_multiplicity": int(peak_basis5.shape[1]),
              "tensor_top_multiplicity": int(peak_basis10.shape[1]),
              "compressed_first_derivative_norm":
                  float(np.linalg.norm(compression, 2)),
              "rows": rows}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result), flush=True)


if __name__ == "__main__":
    main()
