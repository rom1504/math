#!/usr/bin/env python3
"""Replay the joint-field BH mechanism audit; no external data needed.

Run from the repository with .venv/bin/python.  Decimal optimization is a
diagnostic, not an interval certificate.  The .398942... analytic ceiling
and the Paley plateau are proved in the accompanying dated artifact.
"""

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution
from scipy.special import erf


def cube(n):
    bits = np.arange(1 << n, dtype=np.uint64)[:, None]
    return (1 - 2 * ((bits >> np.arange(n, dtype=np.uint64)) & 1).astype(np.int8)).astype(float)


def quadratic(a, x):
    return np.einsum("bi,ij,bj->b", x, a, x, optimize=True) / 2


def paley(p):
    if p % 4 != 1 or any(p % j == 0 for j in range(2, int(p**0.5) + 1)):
        raise ValueError("Use a prime p congruent to 1 modulo 4")
    chi = np.zeros(p, dtype=np.int64)
    for j in range(1, p):
        val = pow(j, (p - 1) // 2, p)
        chi[j] = 1 if val == 1 else -1
    return chi[(np.arange(p)[:, None] - np.arange(p)[None, :]) % p]


def check_paley():
    reports = []
    for p in (5, 13, 17, 29, 37, 101):
        a = paley(p)
        assert np.array_equal(a, a.T)
        assert np.all(np.diag(a) == 0)
        assert np.all(np.abs(a[~np.eye(p, dtype=bool)]) == 1)
        assert np.all(a.sum(axis=0) == 0)
        assert np.array_equal(a @ a, p * np.eye(p, dtype=np.int64) - np.ones((p, p), dtype=np.int64))
        row = {"order": p, "exact_matrix_identities": True}
        if p <= 17:
            x = cube(p)
            q = quadratic(a, x)
            magnetization = x.sum(axis=1)
            envelope = p**0.5 * (p - magnetization**2 / p) / 2
            assert np.max(np.abs(q) - envelope) <= 1e-10
            curve = []
            for t in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2):
                val = float(np.max(np.abs(q + t * p**0.5 * magnetization)) / p**1.5)
                upper = (1 + t * t) / 2 if t <= 1 else t
                assert val <= upper + 1e-12
                if t >= 1:
                    assert abs(val - t) < 1e-12
                curve.append({"field_t": t, "normalized_norm": val, "spectral_upper": upper})
            row["exhaustive_curve"] = curve
        reports.append(row)
    return reports


def check_square_paley_parabola():
    reports = []
    for p in (3, 5, 7, 11, 17):
        chi = np.zeros(p, dtype=np.int64)
        for j in range(1, p):
            chi[j] = 1 if pow(j, (p - 1) // 2, p) == 1 else -1
        nonsquare = int(np.flatnonzero(chi == -1)[0])
        aa = np.tile(np.arange(p), p)
        bb = np.repeat(np.arange(p), p)
        norm = ((aa[:, None] - aa[None, :])**2
                - nonsquare * (bb[:, None] - bb[None, :])**2) % p
        a = chi[norm]
        n = p * p
        assert np.array_equal(a @ a, n * np.eye(n, dtype=np.int64) - np.ones((n, n), dtype=np.int64))
        assert np.all(a.sum(axis=1) == 0)
        row_witnesses = []
        for k in range(p + 1):
            x = np.where(bb < k, -1, 1)
            u = float(x.mean())
            energy = float(x @ a @ x / (2 * n**1.5))
            assert abs(energy - (1 - u * u) / 2) < 1e-12
            row_witnesses.append((u, energy))
        curves = []
        for t in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2):
            witness = max(abs(energy + t * u) for u, energy in row_witnesses)
            upper = (1 + t * t) / 2 if t <= 1 else t
            if t <= 1:
                assert upper - 1 / (2 * p * p) <= witness + 1e-12
            else:
                assert abs(witness - upper) < 1e-12
            curves.append({"t": t, "row_witness_lower": witness, "spectral_upper": upper})
        reports.append({"prime": p, "order": n, "nonsquare": nonsquare,
                        "exact_matrix_and_row_witness_identities": True, "curves": curves})
    return reports


def check_field_lemmas(seed=20260916):
    rng = np.random.default_rng(seed)
    min_hard_slack = float("inf")
    min_soft_slack = float("inf")
    hard_tests = soft_tests = 0
    for n in (3, 4, 5, 6, 7):
        x = cube(n)
        for _ in range(8):
            a = np.triu(rng.choice((-1, 1), (n, n)), 1)
            a = a + a.T
            h = rng.normal(size=n) * n**0.5
            c = float(rng.normal() * n)
            norm = float(np.max(np.abs(c + quadratic(a, x) + x @ h)))
            h_abs = np.abs(h)
            h_total = float(h_abs.sum())
            switch = np.where(h >= 0, 1, -1)
            switched = a * switch[:, None] * switch[None, :]
            for bits in itertools.product((0, 1), repeat=n):
                r = np.asarray(bits, dtype=float)
                vals = quadratic(switched * r[:, None] * r[None, :], x)
                width = (float(vals.max()) - float(vals.min())) / 2
                rhs = h_total - float(h_abs @ r) + width
                slack = norm - rhs
                assert slack >= -1e-10
                min_hard_slack = min(min_hard_slack, slack)
                hard_tests += 1
            for _ in range(40):
                r = rng.uniform(size=n)
                vals = quadratic(switched * r[:, None] * r[None, :], x)
                width = (float(vals.max()) - float(vals.min())) / 2
                rhs = float(h_abs @ (1 - r)) + width
                slack = norm - rhs
                assert slack >= -1e-10
                min_soft_slack = min(min_soft_slack, slack)
                soft_tests += 1
    return {"hard_tests": hard_tests, "soft_tests": soft_tests,
            "minimum_hard_slack": min_hard_slack, "minimum_soft_slack": min_soft_slack}


def check_coordinate_convexity(seed=20260917):
    rng = np.random.default_rng(seed)
    worst_excess = -float("inf")
    for p in (1, 4 / 3, 2):
        for n in (4, 6, 8):
            h = np.abs(rng.normal(size=n))
            weight = rng.uniform(0.1, 2.0)

            def objective(r):
                edge = np.outer(r, r)[np.triu_indices(n, 1)]
                return float(h @ (1 - r) + weight * np.sum(edge**p)**(1 / p))

            vertex_max = max(objective(np.asarray(bits, dtype=float))
                             for bits in itertools.product((0, 1), repeat=n))
            for _ in range(100):
                r = rng.uniform(size=n)
                worst_excess = max(worst_excess, objective(r) - vertex_max)
                assert objective(r) <= vertex_max + 1e-12
                i = int(rng.integers(n))
                lo, hi = r.copy(), r.copy()
                lo[i], hi[i] = 0, 1
                interpolation = (1 - r[i]) * objective(lo) + r[i] * objective(hi)
                assert objective(r) <= interpolation + 1e-12
    return {"p_values": [1, 4 / 3, 2], "largest_sample_interior_excess": worst_excess}


def bootstrap_diagnostics():
    k = float(np.sqrt(2 / np.pi))

    def terms(z):
        alpha, t = z
        p = erf(t / np.sqrt(2))
        linear = k * alpha * np.sqrt(1 - alpha) * np.exp(-t * t / 2)
        residual = (alpha * p)**1.5
        return linear, residual

    def ratio(z):
        linear, residual = terms(z)
        return linear / (1 - residual)

    optimum = differential_evolution(lambda z: -ratio(z), [(1e-7, 1 - 1e-7), (0, 6)],
                                     seed=16, tol=1e-12, polish=True)
    w = 0.4333221116640807
    updated = differential_evolution(lambda z: -(terms(z)[0] + w * terms(z)[1]),
                                    [(1e-7, 1 - 1e-7), (0, 8)], seed=17,
                                    tol=1e-12, polish=True)
    rng = np.random.default_rng(18)
    ceiling_slack = float("inf")
    for _ in range(20000):
        a, p = rng.uniform(size=2)
        lhs = 2 * a * np.sqrt((1 - a) * (1 - p * p)) + (a * p)**1.5
        majorant = np.sqrt(a) * (4 - a) / 3
        assert lhs <= majorant + 1e-12
        assert majorant <= 1 + 1e-12
        ceiling_slack = min(ceiling_slack, 1 - lhs)
    return {"numeric_bootstrap_critical_value": float(-optimum.fun),
            "numeric_critical_alpha_t": optimum.x.tolist(),
            "rigorous_analytic_critical_upper": k / 2,
            "current_half_range_constant": w,
            "numeric_updated_value_with_boundary_cutoff": float(-updated.fun),
            "numeric_updated_alpha_t": updated.x.tolist(),
            "sample_minimum_ceiling_slack": ceiling_slack,
            "qualification": "Optimizer decimals are diagnostics. The analytic ceiling has a written proof."}


def ellipse_bootstrap_checks(seed=2026091611):
    s = float(2 / np.pi)
    critical_a = 1 - s / 3
    critical_b = s / 3 * (1 - s / 3)
    rng = np.random.default_rng(seed)
    max_excess = -float("inf")
    for _ in range(20000):
        baseline = rng.uniform(critical_b, 0.25)
        a = rng.uniform()
        paley_b = 0.25 if a <= 0.5 else a * (1 - a)
        b = min(baseline, paley_b)
        k = a * s
        if k <= 3 * b:
            value = b
        else:
            value = 4 * k**3 / (27 * (k - b)**2)
        max_excess = max(max_excess, value - baseline)
        assert value <= baseline + 1e-12
    return {"necessary_paley_coefficients": "0<=a<=1; b<=1/4 for a<=1/2, b<=a-a^2 for a>=1/2",
            "critical_a": critical_a, "critical_squared_baseline": critical_b,
            "critical_baseline": float(np.sqrt(critical_b)), "sample_checks": 20000,
            "maximum_sample_excess": max_excess,
            "qualification": "The barrier has an analytic proof; random tests are regression only."}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {"status": "PASS", "paley": check_paley(),
              "square_paley_parabola": check_square_paley_parabola(),
              "field_lemmas": check_field_lemmas(),
              "coefficient_soft_reserve": check_coordinate_convexity(),
              "bootstrap": bootstrap_diagnostics(),
              "ellipse_bootstrap": ellipse_bootstrap_checks()}
    rendered = json.dumps(report, indent=2, sort_keys=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n")
    print(rendered)


if __name__ == "__main__":
    main()
