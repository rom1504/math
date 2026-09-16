#!/usr/bin/env python3
"""Exact signing-filter regression and numerical generic moment matching.

The signing tests use integer/Fraction arithmetic.  Quadrature tests are
floating-point regression, not interval certificates; the accompanying
artifact proves their positivity and exactness algebraically.
"""

import argparse
import json
from fractions import Fraction
from math import acosh, comb, cos, log, pi
from pathlib import Path

import numpy as np
from scipy.optimize import brentq


def signing_values(n, rng):
    a = np.triu(rng.choice((-1, 1), (n, n)), 1)
    a = a + a.T
    x = 1 - 2 * ((np.arange(1 << n)[:, None] >> np.arange(n)) & 1)
    return [int(v) for v in np.einsum("bi,ij,bj->b", x, a, x) // 2]


def integer_filter_coefficients(degree, cap):
    """Coefficients of cap**degree * T_degree(X/cap)."""
    old, current = [1], [0, 1]
    if degree == 0:
        return old
    for _ in range(1, degree):
        new = [0] + [2 * v for v in current]
        for j, v in enumerate(old):
            new[j] -= cap * cap * v
        old, current = current, new
    return current


def evaluate_polynomial(coefficients, x):
    value = 0
    for a in reversed(coefficients):
        value = value * x + a
    return value


def exact_signing_checks(seed=2026091621):
    rng = np.random.default_rng(seed)
    reports = []
    for n in (5, 6, 8, 9):
        energies = signing_values(n, rng)
        true_cap = max(abs(v) for v in energies)
        endpoint_count = sum(abs(v) == true_cap for v in energies)
        endpoint_mass = Fraction(endpoint_count, 1 << n)
        assert endpoint_mass >= Fraction(1, 1 << (n - 1))
        moments = [Fraction(sum(v**j for v in energies), 1 << n) for j in range(33)]
        assert moments[2] == Fraction(n * (n - 1), 2)
        rows = []
        for cap in sorted({max(1, true_cap * 4 // 5), true_cap}):
            for degree in (1, 2, 4, 8, 16):
                coefficients = integer_filter_coefficients(degree, cap)
                values = [evaluate_polynomial(coefficients, v) for v in energies]
                direct = Fraction(sum(v * v for v in values), (1 << n) * cap**(2 * degree))
                square = [0] * (2 * degree + 1)
                for i, a in enumerate(coefficients):
                    for j, b in enumerate(coefficients):
                        square[i + j] += a * b
                from_moments = sum(a * moments[j] for j, a in enumerate(square)) / cap**(2 * degree)
                assert direct == from_moments
                endpoint_value = evaluate_polynomial(coefficients, true_cap)
                endpoint_lower = endpoint_mass * Fraction(endpoint_value**2, cap**(2 * degree))
                assert direct >= endpoint_lower
                if cap >= true_cap:
                    assert direct <= 1
                guaranteed_detect = Fraction(endpoint_value**2, (1 << (n - 1)) * cap**(2 * degree)) > 1
                if guaranteed_detect:
                    assert direct > 1
                rows.append({"candidate_cap": cap, "degree": degree,
                             "exact_L2_squared": str(direct), "L2_detects": direct > 1,
                             "universal_endpoint_guarantee_detects": guaranteed_detect})
        reports.append({"n": n, "true_cap": true_cap, "endpoint_mass": str(endpoint_mass),
                        "tests": rows})
    return reports


def chebyshev_values(x, degree):
    x = np.asarray(x, dtype=float)
    values = [np.ones_like(x)]
    if degree:
        values.append(x.copy())
    for _ in range(1, degree):
        values.append(2 * x * values[-1] - values[-2])
    return np.asarray(values)


def generic_quadrature_checks():
    reports = []
    for degree, atom_exponent in ((4, 20), (8, 30), (16, 40), (32, 40), (64, 40)):
        p = 2.0**(-atom_exponent)
        ell = -log(p)
        assert p < 1 / (2 * degree + 1)

        def kernel_at_alpha(alpha):
            return 1 + 2 * np.sum(np.cosh(np.arange(1, degree + 1) * alpha)**2)

        alpha = brentq(lambda z: log(kernel_at_alpha(z)) - ell,
                       0, (ell + log(2)) / (2 * degree), xtol=1e-14)
        radius = float(np.cosh(alpha))
        r_values = chebyshev_values(radius, degree)
        coefficients = 2 * r_values
        coefficients[0] = 1
        roots = np.polynomial.chebyshev.chebroots(coefficients)
        if np.iscomplexobj(roots):
            assert np.max(np.abs(roots.imag)) < 1e-9
            roots = roots.real
        assert np.all(np.abs(roots) < 1)
        nodes = np.append(roots, radius)
        values = chebyshev_values(nodes, degree)
        diagonal = 1 + 2 * np.sum(values[1:]**2, axis=0)
        weights = 1 / diagonal
        assert np.all(weights > 0)
        assert np.min(weights[:-1]) >= 1 / (2 * degree + 1) - 1e-12
        assert abs(weights[-1] / p - 1) < 1e-9
        base_nodes = np.asarray([cos((2 * j - 1) * pi / (2 * (degree + 1)))
                                 for j in range(1, degree + 2)])
        max_error = 0.0
        max_between_error = 0.0
        for k in range(2 * degree + 1):
            expected = 0.0 if k % 2 else comb(k, k // 2) / 2.0**k
            computed = float(weights @ nodes**k)
            base = float(np.mean(base_nodes**k))
            max_error = max(max_error, abs(computed - expected))
            max_between_error = max(max_between_error, abs(computed - base))
        assert max_error < 1e-8
        assert max_between_error < 1e-8
        lower_alpha = (ell - log(2 * degree + 1)) / (2 * degree)
        upper_alpha = (ell + log(2)) / (2 * degree)
        assert lower_alpha <= alpha + 1e-12 <= upper_alpha + 1e-12
        reports.append({"degree": degree, "minimum_atom_target": p, "outer_radius": radius,
                        "base_cap": float(np.max(np.abs(base_nodes))), "alpha": alpha,
                        "alpha_lower": lower_alpha, "alpha_upper": upper_alpha,
                        "outer_atom_weight": float(weights[-1]),
                        "minimum_interior_weight": float(np.min(weights[:-1])),
                        "sum_weights": float(weights.sum()), "maximum_moment_error": max_error,
                        "maximum_two_measure_moment_difference": max_between_error,
                        "qualification": "Floating-point generic-measure regression; NOT an actual signing or a dyadic histogram."})
    return reports


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = {"status": "PASS", "exact_actual_signing_filters": exact_signing_checks(),
              "generic_moment_matching": generic_quadrature_checks(),
              "archive_status": "Chebyshev upper recovery already proved in blank_slate_direct_attack_2026_08_21.md Section 2.1."}
    rendered = json.dumps(report, indent=2, sort_keys=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n")
    print(rendered)


if __name__ == "__main__":
    main()
