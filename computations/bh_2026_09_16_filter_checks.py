"""Exact Walsh transforms of nonlinear quadratic-energy filters.

All matrices, energies, transforms, and Chebyshev identities are integer
certificates. Coefficient norms use mpmath and are explicitly diagnostics,
not interval certificates. The asymptotic claims are proved in the artifacts.
"""

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def fwht(values):
    out = list(map(int, values))
    width = 1
    while width < len(out):
        for start in range(0, len(out), 2 * width):
            for j in range(start, start + width):
                a, b = out[j], out[j + width]
                out[j], out[j + width] = a + b, a - b
        width *= 2
    return out


def regular_hadamard(k):
    base = np.ones((4, 4), dtype=np.int64) - 2 * np.eye(4, dtype=np.int64)
    h = np.ones((1, 1), dtype=np.int64)
    for _ in range(k):
        h = np.kron(h, base)
    n = len(h)
    assert np.array_equal(h @ h, n * np.eye(n, dtype=np.int64))
    delta = (-1) ** k
    assert np.all(np.diag(h) == delta)
    return h - delta * np.eye(n, dtype=np.int64)


def energies(a):
    n = len(a)
    spins = 1 - 2 * ((np.arange(2**n, dtype=np.int64)[:, None]
                      >> np.arange(n, dtype=np.int64)) & 1)
    doubled = np.sum((spins @ a) * spins, axis=1)
    assert np.all(doubled % 2 == 0)
    return list(map(int, doubled // 2))


def mp_norm(histogram, exponent):
    return mp.fsum(mp.mpf(v)**exponent * count
                   for v, count in histogram.items())**(1 / exponent)


def check_case(name, a, degrees, regular=False):
    n = len(a)
    root = math.isqrt(n)
    assert root * root == n
    assert np.array_equal(a, a.T) and np.all(np.diag(a) == 0)
    assert set(map(int, a[np.triu_indices(n, 1)])) <= {-1, 1}
    energy = energies(a)
    cap = max(map(abs, energy))
    if regular:
        assert cap == n * (root + 1) // 2
    # c=2/5.  Y/c = (5 H)/(2 n sqrt(n)).
    denominator = 2 * n * root
    f_previous, f_current = [1] * len(energy), [5 * h for h in energy]
    reports = []
    for d in range(1, max(degrees) + 1):
        if d in degrees:
            scale = denominator**d
            transform = fwht(f_current)
            assert all(v % len(energy) == 0 for v in transform)
            coefficients = [v // len(energy) for v in transform]
            assert sum(v*v for v in f_current) == len(energy) * sum(
                v*v for v in coefficients)
            assert fwht(coefficients) == f_current
            levels = {}
            for mask, coefficient in enumerate(coefficients):
                if coefficient:
                    r = bin(mask).count("1")
                    assert r % 2 == 0 and r <= 2*d
                    levels.setdefault(r, Counter())[abs(coefficient)] += 1
            if 2*d <= n:
                assert max(levels) == 2*d
            p = mp.mpf(4*d) / (2*d+1)
            global_norm = mp_norm(Counter(map(abs, coefficients)), p) / scale
            weighted = {}
            for s in [mp.mpf(5), mp.mpf(5)/2]:
                squares = []
                for r, histogram in levels.items():
                    if r:
                        q = mp.mpf(2*r) / (r+1)
                        squares.append(mp_norm(histogram, q)**2 / r**(2*s))
                weighted[str(s)] = mp.nstr(mp.sqrt(mp.fsum(squares))/scale, 30)
            reports.append({
                "degree": d, "actual_Walsh_degree": max(levels),
                "denominator": str(scale),
                "exact_max_abs_numerator": str(max(map(abs, f_current))),
                "exact_nonzero_coefficients": sum(sum(h.values()) for h in levels.values()),
                "exact_level_histograms_abs_numerators": {
                    str(r): {str(v): count for v, count in sorted(hist.items())}
                    for r, hist in sorted(levels.items())},
                "numeric_coefficient_norm": mp.nstr(global_norm, 30),
                "numeric_weighted_norms": weighted,
                "numeric_function_supremum": mp.nstr(
                    mp.mpf(max(map(abs, f_current)))/scale, 30),
            })
        f_previous, f_current = f_current, [
            10*h*now - denominator**2*before
            for h, now, before in zip(energy, f_current, f_previous)]
    return {"name": name, "n": n, "matrix": a.tolist(),
            "matrix_sha256": hashlib.sha256(a.tobytes()).hexdigest(),
            "cap": cap, "normalized_cap": mp.nstr(mp.mpf(cap)/(n*root), 30),
            "exact_energy_histogram": dict(sorted(Counter(energy).items())),
            "filters": reports}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    mp.mp.dps = 60
    witness = json.loads((ROOT / "computations/results/exact_m9.json").read_text())
    cases = [
        check_case("regular_Hadamard_4", regular_hadamard(1), [1, 2, 3, 4, 6], True),
        check_case("regular_Hadamard_16", regular_hadamard(2), [1, 2, 3, 4, 6, 8, 12], True),
        check_case("stored_cap12_order9", np.array(witness["matrix"], dtype=np.int64),
                   [1, 2, 3, 4, 6]),
        check_case("all_positive_16", np.ones((16, 16), dtype=np.int64)
                   - np.eye(16, dtype=np.int64), [1, 2, 3, 4, 6, 8, 12]),
    ]
    result = {"status": "PASS", "normalization": "P(Y)=T_d(Y/(2/5)), Y=H/n^(3/2)",
              "evidence": "Integer certificates plus separately labeled numerical norms; not asymptotic evidence",
              "source_dependency": "computations/results/exact_m9.json (cap replayed; optimality unused)",
              "cases": cases}
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"status": "PASS", "cases": len(cases),
                      "filters": sum(len(c["filters"]) for c in cases),
                      "caps": {c["name"]: c["cap"] for c in cases}}))


if __name__ == "__main__":
    main()
