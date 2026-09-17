"""Checks for the exact-new-spin composition criteria (artifact §§11–12).

Run .venv/bin/python computations/paper_bernoulli_2026_09_17_abs_columns.py.
Numerics corroborate the displayed proofs; no random-matrix claim is inferred.
"""

import json
import math

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar


def rate(w):
    return (w * w - 1.0 - 2.0 * math.log(w)) / 2.0


def chi_threshold(k_value, inner):
    r = math.sqrt(inner)
    if k_value == 0.0:
        return r
    k = k_value / r
    w = (k + math.sqrt(k * k + 4.0)) / 2.0
    return r * (w / 2.0 + math.asinh(k / 2.0) / k)


def union_threshold(k_value, inner):
    ell = math.sqrt(2.0 * math.log(2.0) * inner)
    if k_value <= ell:
        return ell
    return (k_value * k_value + ell * ell) / (2.0 * k_value)


def mixed_resource(k_value, inner, tau):
    return (2.0 * k_value / tau + inner / tau ** 2
            + math.pi / 2.0 * math.sqrt(k_value * inner) / tau ** 1.5)


def discrete_resource(epsilon, inner, k_value, alpha, tau, ratio=1.015):
    eta0 = 0.01
    levels = [eta0]
    while levels[-1] > epsilon ** 2:
        levels.append(levels[-1] / ratio)
    eta = np.array(levels)
    g = inner + 0.4 * eta ** alpha
    s = k_value * eta / g
    assert np.all(g < 1.0)
    assert np.all(s < math.log(2.0))
    assert np.all(np.diff(g) <= 0.0)
    assert np.all(np.diff(s) <= 0.0)
    denominators = (eta / ratio + tau * epsilon) ** 2
    denominators[-1] = (tau * epsilon) ** 2
    precision = (epsilon ** 2 + 2.0 * epsilon ** 1.5 * np.sqrt(s)
                 + 2.0 * epsilon * s) / denominators
    running = np.maximum.accumulate(precision)
    increments = np.diff(np.concatenate(([0.0], running)))
    return float(g @ increments)


def main():
    # Explicit square-MGF consequence: J(sqrt(1+2sqrt(u)+2u)) >= u.
    tail_checks = 0
    min_tail_gap = float("inf")
    for u in np.logspace(-10, 7, 600):
        w = math.sqrt(1.0 + 2.0 * math.sqrt(u) + 2.0 * u)
        gap = rate(w) - u
        assert gap >= -2e-13 * max(1.0, u)
        min_tail_gap = min(min_tail_gap, gap)
        tail_checks += 1

    beta_integral, error = quad(lambda u: math.sqrt(u) / (1.0 + u) ** 3,
                                0.0, np.inf, epsabs=1e-11)
    assert abs(beta_integral - math.pi / 8.0) < 1e-10

    # Independently maximize the scalar variational problems.
    threshold_checks = []
    for inner in (0.04, 0.2, 0.7):
        for k_value in (0.01, 0.1, 0.4, 1.2):
            r = math.sqrt(inner)
            k = k_value / r
            optimum = minimize_scalar(lambda w: -(w - rate(w) / k),
                                      bounds=(1.0, max(5.0, 3.0 * k + 3.0)),
                                      method="bounded",
                                      options={"xatol": 1e-12})
            formula = chi_threshold(k_value, inner)
            numerical = -r * optimum.fun
            assert abs(formula - numerical) < 1e-10
            tau = formula * 1.001
            for z in np.concatenate(([0.0], np.logspace(-8, 5, 300))):
                assert rate(tau / r + z) > k * z

            ell = math.sqrt(2.0 * math.log(2.0) * inner)
            union_opt = minimize_scalar(
                lambda u: -(math.sqrt(2.0 * k_value * u + ell ** 2) - u),
                bounds=(0.0, max(5.0, 3.0 * k_value + 3.0)),
                method="bounded", options={"xatol": 1e-12})
            numerical_union = max(ell, -union_opt.fun)
            exact_union = union_threshold(k_value, inner)
            assert abs(numerical_union - exact_union) < 1e-10
            threshold_checks.append({
                "e0": inner, "K": k_value,
                "F_chi": formula, "F_union": exact_union,
                "best": min(formula, exact_union),
            })

    budgets = []
    for inner, k_value, alpha, tau in ((0.0, 0.2, 0.5, 0.7),
                                     (0.15, 0.04, 0.4, 0.7),
                                     (0.3, 0.2, 0.25, 1.1)):
        values = [{"epsilon": eps,
                   "resource": discrete_resource(eps, inner, k_value, alpha, tau)}
                  for eps in (1e-3, 1e-5, 1e-7, 1e-9)]
        # The theorem is an asymptotic upper bound, not a finite-epsilon one.
        bound = mixed_resource(k_value, inner, tau)
        assert values[-1]["resource"] < bound + 0.02
        budgets.append({"e0": inner, "K": k_value, "alpha": alpha,
                        "tau": tau, "asymptotic_upper": bound,
                        "finite_resources": values})
    print(json.dumps({"status": "PASS", "tail_checks": tail_checks,
                      "minimum_tail_gap": min_tail_gap,
                      "beta_integral": beta_integral,
                      "threshold_checks": threshold_checks,
                      "resource_checks": budgets}, indent=2))


if __name__ == "__main__":
    main()
