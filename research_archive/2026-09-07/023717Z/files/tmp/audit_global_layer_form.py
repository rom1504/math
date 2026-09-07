#!/usr/bin/env python3
"""Audit homogeneous Fourier-layer compressions of the signing Hamiltonian."""

from __future__ import annotations

import itertools
import json
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import eigsh


ROOT = Path("/home/math/quadra")


def cap(a: np.ndarray) -> int:
    n = len(a)
    best = 0
    for code in range(1 << (n - 1)):
        x = np.ones(n, dtype=np.int64)
        for i in range(1, n):
            if code & (1 << (i - 1)):
                x[i] = -1
        best = max(best, abs(int(x @ a @ x) // 2))
    return best


def layer_operator(a: np.ndarray, k: int) -> sparse.csr_matrix:
    n = len(a)
    subsets = list(itertools.combinations(range(n), k))
    index = {s: i for i, s in enumerate(subsets)}
    rows: list[int] = []
    cols: list[int] = []
    vals: list[int] = []
    for row, s_tuple in enumerate(subsets):
        s = set(s_tuple)
        for i in s:
            for j in range(n):
                if j not in s:
                    t = tuple(sorted((s - {i}) | {j}))
                    rows.append(row)
                    cols.append(index[t])
                    vals.append(int(a[i, j]))
    return sparse.csr_matrix((vals, (rows, cols)), shape=(len(subsets), len(subsets)))


def self_complementary(n: int, seed: int) -> np.ndarray:
    assert n % 4 == 0
    r = n // 4
    rng = np.random.default_rng(seed)
    c0 = np.zeros((r, r), dtype=np.int8)
    c2 = np.zeros((r, r), dtype=np.int8)
    for c, diagonal in ((c0, False), (c2, True)):
        for i in range(r):
            if diagonal:
                c[i, i] = rng.choice([-1, 1])
            for j in range(i + 1, r):
                c[i, j] = c[j, i] = rng.choice([-1, 1])
    c1 = rng.choice([-1, 1], size=(r, r)).astype(np.int8)
    blocks = [c0, c1, c2, -c1.T]
    a = np.zeros((n, n), dtype=np.int8)
    for layer in range(4):
        for d in range(4):
            target = (layer + d) % 4
            a[layer*r:(layer+1)*r, target*r:(target+1)*r] = ((-1) ** layer) * blocks[d]
    np.fill_diagonal(a, 0)
    assert np.array_equal(a, a.T)
    return a


def load_matrix(path: str, key: str = "matrix") -> np.ndarray:
    data = json.loads((ROOT / path).read_text())
    return np.asarray(data[key], dtype=np.int8)


def audit(name: str, a: np.ndarray) -> dict[str, object]:
    n = len(a)
    k = n // 2
    t = layer_operator(a, k)
    d = t.shape[0]
    degree = k * (n-k)
    assert t.nnz == d * degree
    assert (t - t.T).nnz == 0
    lam = eigsh(t.astype(float), k=2, which="BE", return_eigenvectors=False, tol=1e-11)
    op = float(max(abs(lam)))
    t2 = t @ t
    tr2 = int(np.dot(t.data.astype(np.int64), t.data.astype(np.int64)))
    tr4 = int(np.dot(t2.data.astype(np.int64), t2.data.astype(np.int64)))
    m = cap(a)
    coeff = 2*k*(n-k)/(n*(n-1))
    return {
        "name": name,
        "n": n,
        "k": k,
        "dimension": d,
        "degree": degree,
        "cap": m,
        "layer_operator_norm": op,
        "rayleigh_guarantee": coeff*m,
        "operator_norm_over_cap": op/m,
        "second_moment_certificate": (tr2/d)**0.5,
        "fourth_moment_certificate": (tr4/d)**0.25,
        "trace2": tr2,
        "trace4": tr4,
        "trace4_over_dimension_n4": tr4/(d*n**4),
        "aligned_quarter_turn_self_complementarity_error": (
            int(np.max(np.abs(np.roll(np.roll(a, -n//4, axis=0), -n//4, axis=1)+a)))
            if name.startswith("selfcomp_") else None
        ),
        "trace_A4_over_n3": float(
            np.trace(
                (a.astype(np.int64) @ a.astype(np.int64))
                @ (a.astype(np.int64) @ a.astype(np.int64))
            )
        )/n**3,
    }


def main() -> None:
    cases = [
        ("exact_m6_conference", load_matrix("computations/results/exact_m6.json")),
        ("exact_m8_class0", np.asarray(json.loads((ROOT/"computations/results/m8_minimizer_orbits.json").read_text())["classes"][0]["representative_matrix"], dtype=np.int8)),
        ("exact_m8_class1_nonconference", np.asarray(json.loads((ROOT/"computations/results/m8_minimizer_orbits.json").read_text())["classes"][1]["representative_matrix"], dtype=np.int8)),
        ("exact_m10", load_matrix("computations/results/exact_m10.json")),
        ("conference_gf9_n10", load_matrix("computations/results/conference_order10_gf9.json", "conference_matrix")),
        ("heuristic_m12", load_matrix("computations/results/heuristic_m12_seed20260731.json")),
        ("exact_m14_conference", load_matrix("computations/results/conference_completion_m13.json", "conference_matrix")),
        ("selfcomp_wigner_n8_seed20260813", self_complementary(8, 20260813)),
        ("selfcomp_wigner_n12_seed20260813", self_complementary(12, 20260813)),
    ]
    results=[]
    for name,a in cases:
        row=audit(name,a)
        results.append(row)
        print(json.dumps(row,sort_keys=True),flush=True)
    print(json.dumps({"schema":"global-layer-form-audit-v1","results":results},indent=2,sort_keys=True))


if __name__ == "__main__":
    main()
