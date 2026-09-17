"""Audit computations for cumulative-precision Gram--Schmidt augmentation.

Run from the repository root:
    .venv/bin/python computations/paper_discrepancy_2026_09_17_precision.py

These numerical checks corroborate the proof in the matching artifact;
they do not replace the all-order mathematical argument.
"""

from __future__ import annotations

import json
import math

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq


def entropy(g: float) -> float:
    if g == 0.0 or g == 1.0:
        return 0.0
    return -g * math.log(g) - (1.0 - g) * math.log1p(-g)


def complexity(eta: float, kappa: float = 1.0) -> float:
    """Solve g h(g) = kappa eta, using a log-coordinate near zero."""
    if eta == 0.0:
        return 0.0
    if kappa * eta >= 0.5 * math.log(2.0):
        raise ValueError("eta is outside the small-complexity regime")
    log_target = math.log(kappa * eta)

    def residual(log_g: float) -> float:
        g = math.exp(log_g)
        return log_g + math.log(entropy(g)) - log_target

    log_g = brentq(residual, min(-2.0, log_target - 3.0), math.log(0.5))
    return math.exp(log_g)


def rho(g: float) -> float:
    if g == 0.0:
        return 0.5
    derivative = math.log1p(-g) - math.log(g)
    return g * derivative / (entropy(g) + g * derivative)


def principal_resource(
    epsilon: float, *, tau: float, ratio: float, eta0: float, kappa: float = 1.0
) -> float:
    """Exact integral (13), transformed by u=t/(1-t) for stability."""
    cap = tau * epsilon
    tmax = eta0 / (eta0 + ratio * cap)

    def integrand(t: float) -> float:
        eta = ratio * cap * t / (1.0 - t)
        return max(0.0, 2.0 * t - rho(complexity(eta, kappa)))

    val, error = quad(integrand, 0.0, tmax, epsabs=2e-10, limit=250)
    assert error < 2e-7
    return 2.0 * kappa * ratio / tau * val


def full_integral_bound(
    epsilon: float, *, tau: float, ratio: float, eta0: float, kappa: float = 1.0
) -> float:
    """Integral in (12), with the new-word entropy term retained."""
    cap = tau * epsilon
    tmax = eta0 / (eta0 + ratio * cap)

    def integrand(t: float) -> float:
        eta = ratio * cap * t / (1.0 - t)
        g = complexity(eta, kappa)
        main = 2.0 * kappa * ratio / tau * (2.0 * t - rho(g))
        extra = 4.0 * math.log(2.0) / (tau * tau) * g * (1.0 - t)
        return max(0.0, main + extra)

    val, error = quad(integrand, 0.0, tmax, epsabs=2e-10, limit=250)
    assert error < 2e-7
    return val


def discrete_resource(
    epsilon: float, *, tau: float, ratio: float, eta0: float, kappa: float = 1.0
) -> dict[str, float | int]:
    levels = [eta0]
    while levels[-1] > epsilon * epsilon:
        levels.append(levels[-1] / ratio)
    g = np.array([complexity(e, kappa) for e in levels])
    b = np.array(
        [
            2.0 * epsilon * (entropy(gg) + epsilon * math.log(2.0))
            / ((e / ratio if j + 1 < len(levels) else 0.0) + tau * epsilon) ** 2
            for j, (e, gg) in enumerate(zip(levels, g))
        ]
    )
    running = np.maximum.accumulate(b)
    increments = np.diff(np.concatenate(([0.0], running)))
    resource = float(g @ increments)
    exclusive = float(g @ b)
    boundary = float(g[0] * b[0])
    terminal = float(g[-1] * b[-1])
    upper_bound = boundary + terminal + full_integral_bound(
        epsilon, tau=tau, ratio=ratio, eta0=eta0, kappa=kappa
    )
    assert resource <= upper_bound + 3e-8
    assert resource <= exclusive + 1e-10
    return {
        "epsilon": epsilon,
        "levels": len(levels),
        "cumulative_resource": resource,
        "exclusive_resource": exclusive,
        "integral_upper_bound": upper_bound,
    }


def audit_bellman(seed: int = 910172026, trials: int = 5000) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    worst = 0.0
    for _ in range(trials):
        x = float(rng.uniform(-0.999999, 0.999999))
        t = float(rng.uniform(-20.0, 20.0))
        # Derivative of centered-sign log MGF, using a stable hyperbolic form.
        psi_prime = math.tanh(t + math.atanh(x)) - x
        gap = t * t - psi_prime * psi_prime
        worst = min(worst, gap)
        assert gap >= -1e-12
        a = abs(t) + float(rng.uniform(0.0, 2.0))
        b = t * t / a + float(rng.uniform(0.0, 2.0))
        determinant = a * b - psi_prime * psi_prime
        assert determinant >= -1e-12
    return {"min_t_squared_minus_derivative_squared": worst}


def audit_harmonic(seed: int = 17092026) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    n = 17
    g = rng.uniform(0.05, 1.0, size=8)
    a_squared = rng.uniform(0.05, 0.9, size=8)
    precision = a_squared / g
    coefficients = precision / precision.sum()
    actual_norm_squared = n * np.sum(coefficients * coefficients * g / a_squared)
    predicted = n / precision.sum()
    assert abs(actual_norm_squared - predicted) < 1e-12
    one_block = n / precision.max()
    return {
        "combined_test_norm_squared": float(actual_norm_squared),
        "harmonic_formula": float(predicted),
        "best_single_block_norm_squared": float(one_block),
    }


def audit_boolean_orthogonality(n: int = 7) -> dict[str, int | float]:
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    words = np.array(
        [[1.0] + [1.0 if mask & (1 << j) else -1.0 for j in range(n - 1)]
         for mask in range(1 << (n - 1))]
    )
    query = np.column_stack([words[:, i] * words[:, j] for i, j in edges])
    error = float(np.max(np.abs(query.T @ query / len(words) - np.eye(len(edges)))))
    assert error == 0.0
    return {"n": n, "edges": len(edges), "orthogonality_error": error}


def audit_bias_obstruction() -> list[dict[str, float]]:
    rows = []
    for p in [1e-1, 1e-2, 1e-4, 1e-8]:
        b = 1.0 - 2.0 * p
        lam = -2.0 * math.log(1.0 / p) / (1.0 + b)
        log_mgf = float(np.logaddexp(
            math.log(p) + lam * (-1.0 - b),
            math.log1p(-p) + lam * (1.0 - b),
        ))
        lower = (1.0 + b) ** 2 / (2.0 * math.log(1.0 / p))
        proxy_at_lambda = 2.0 * log_mgf / (lam * lam)
        variance = 4.0 * p * (1.0 - p)
        assert proxy_at_lambda >= lower * (1.0 - 1e-12)
        rows.append({"p": p, "variance": variance,
                     "proxy_lower_bound": lower, "ratio": lower / variance})
    return rows


def main() -> None:
    exact_integral, _ = quad(lambda t: max(0.0, 2.0 * t - 0.5), 0.0, 1.0)
    assert abs(exact_integral - 9.0 / 16.0) < 1e-12
    parameters = {"tau": 1.4, "ratio": 1.08, "eta0": 0.02}
    principal = []
    for epsilon in [1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-14]:
        value = principal_resource(epsilon, **parameters)
        principal.append({"epsilon": epsilon, "integral": value})
    phi = (1.0 + math.sqrt(5.0)) / 2.0
    report = {
        "status": "all deterministic checks passed",
        "constants": {
            "bellman_integral": exact_integral,
            "new_coefficient": 9.0 / 8.0,
            "archived_coefficient": 20.0 * phi ** 5,
            "only_sharpen_subgaussian": phi ** 5 / 2.0,
            "new_liminf_complexity_factor": 4.0 / 3.0,
        },
        "bellman": audit_bellman(),
        "harmonic": audit_harmonic(),
        "boolean": audit_boolean_orthogonality(),
        "bias": audit_bias_obstruction(),
        "resource_parameters": parameters,
        "principal_integrals": principal,
        "principal_limit": 9.0 * parameters["ratio"] / (8.0 * parameters["tau"]),
        "discrete_resources": [
            discrete_resource(epsilon, **parameters)
            for epsilon in [1e-3, 1e-5, 1e-7, 1e-9]
        ],
    }
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
