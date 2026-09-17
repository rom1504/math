#!/usr/bin/env python3
"""Exact low-order noise polynomials for Q on the coefficient-sign cube.

All Walsh transforms use integer arithmetic. Derivative signs are first scanned
numerically and then certified by exact rational polynomial root isolation for
each distinct minimizer polynomial. This is a diagnostic, not an asymptotic
theorem. No outputs are written outside the repository-local working folder.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import sympy as sp


def fwht(values: np.ndarray) -> np.ndarray:
    out = values.copy()
    half = 1
    while half < len(out):
        blocks = out.reshape(-1, half * 2)
        left = blocks[:, :half].copy()
        right = blocks[:, half:].copy()
        blocks[:, :half] = left + right
        blocks[:, half:] = left - right
        half *= 2
    return out


def calculate(n: int) -> dict:
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    d = len(edges)
    count = 1 << d
    masks = np.arange(count, dtype=np.uint64)
    degree = np.zeros(count, dtype=np.uint8)
    for e in range(d):
        degree += ((masks >> e) & 1).astype(np.uint8)
    vertex_parity = np.zeros(count, dtype=np.uint16)
    for e, (i, j) in enumerate(edges):
        vertex_parity ^= (((masks >> e) & 1) * ((1 << i) | (1 << j))).astype(np.uint16)
    spin_masks = np.arange(1 << (n - 1), dtype=np.uint64)
    spins = np.ones((len(spin_masks), n), dtype=np.int64)
    for i in range(1, n):
        spins[:, i] = 1 - 2 * ((spin_masks >> (i - 1)) & 1).astype(np.int64)
    cuts = np.stack([spins[:, i] * spins[:, j] for i, j in edges])
    values = np.empty(count, dtype=np.int64)
    for start in range(0, count, 8192):
        current = masks[start:start + 8192]
        signs = np.stack([1 - 2 * ((current >> e) & 1).astype(np.int64) for e in range(d)], axis=1)
        values[start:start + len(current)] = np.max(np.abs(signs @ cuts), axis=1)
    transform = fwht(values)
    forbidden = (vertex_parity != 0) | (degree % 2 != 0)
    assert np.all(transform[forbidden] == 0)
    levels = np.zeros((d + 1, count), dtype=np.int64)
    for k in range(d + 1):
        if np.any(transform[degree == k]):
            levels[k] = fwht(np.where(degree == k, transform, 0))
    assert np.array_equal(levels.sum(axis=0), count * values)
    min_value = int(values.min())
    minimizers = np.flatnonzero(values == min_value)
    unique = np.unique(levels[:, minimizers].T, axis=0)
    q = sp.Symbol("q", real=True)
    polynomials = []
    for coefficient in unique:
        # All powers are even: q=t^2. Derivative in q has the same sign as
        # derivative in t for 0<t<1.
        polynomial = sp.Poly(sum(int(coefficient[k]) * q ** (k // 2)
                                 for k in range(0, d + 1, 2)), q)
        derivative = polynomial.diff()
        intervals = sp.polys.polytools.intervals(derivative, eps=sp.Rational(1, 10 ** 10)) if derivative else []
        roots_inside = [(a, b, multiplicity) for (a, b), multiplicity in intervals
                        if b > 0 and a < 1]
        endpoints = [sp.Rational(0), sp.Rational(1)]
        for a, b, _ in roots_inside:
            endpoints.extend((max(sp.Rational(0), a), min(sp.Rational(1), b)))
        endpoints = sorted(set(endpoints))
        test_points = endpoints + [(a + b) / 2 for a, b in zip(endpoints[:-1], endpoints[1:])]
        signs = [int(sp.sign(derivative.eval(point))) for point in test_points]
        polynomials.append({"numerator_coefficients_t": coefficient.tolist(),
                            "denominator": count,
                            "polynomial_q": str(polynomial.as_expr() / count),
                            "derivative_roots_in_unit_interval": [[str(a), str(b), m] for a, b, m in roots_inside],
                            "derivative_nonpositive_exact_root_audit": all(sign <= 0 for sign in signs),
                            "derivative_test_signs": signs})
    c4_coefficients = np.unique(transform[degree == 4])
    c4_nonzero = sorted(set(int(v) for v in transform[(degree == 4) & (vertex_parity == 0)]))
    # Look for failures of monotonicity over the entire signing class, and for
    # curvature signs among minimizers. These scans are explicitly numerical.
    grid = np.linspace(0.001, 0.999, 101)
    first_increasing = None
    for t in grid:
        derivative = np.sum(np.arange(d + 1)[:, None] * levels * np.power(t, np.maximum(np.arange(d + 1) - 1, 0))[:, None], axis=0)
        violating = np.flatnonzero(derivative > 1e-7)
        if len(violating):
            first_increasing = {"edge_mask": int(violating[0]), "cap": int(values[violating[0]]),
                                "t": float(t), "derivative": float(derivative[violating[0]] / count)}
            break
    return {"n": n, "edge_dimension": d, "number_signings": count,
            "minimum_cap": min_value, "minimizer_count": len(minimizers),
            "mean_random_cap_exact": str(sp.Rational(int(transform[0]), count)),
            "nonzero_four_cycle_coefficient_numerators": c4_nonzero,
            "four_cycle_denominator": count,
            "four_cycle_all_degree4_coefficients": c4_coefficients.tolist(),
            "forbidden_fourier_support_exactly_zero": True,
            "minimizer_noise_polynomials": polynomials,
            "first_nonminimizer_increasing_noise_curve_grid_witness": first_increasing}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orders", type=int, nargs="+", default=[3, 4, 5, 6])
    parser.add_argument("--output", type=Path,
                        default=Path("tmp/paper_portfolio_2026_09_17/localization/noise_polynomials.json"))
    args = parser.parse_args()
    results = []
    for n in args.orders:
        result = calculate(n)
        results.append(result)
        print(json.dumps(result, indent=2), flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(results, indent=2) + "\n")


if __name__ == "__main__":
    main()
