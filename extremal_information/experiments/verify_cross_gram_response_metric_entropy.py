#!/usr/bin/env python3
"""Finite diagnostics for cross-Gram response metric entropy."""

from __future__ import annotations

from itertools import product
from math import sqrt

import numpy as np
from scipy.optimize import minimize_scalar


def q_metric(K: np.ndarray, L: np.ndarray) -> float:
    p = len(K)
    best = 0.0
    for eps0 in product((-1.0, 1.0), repeat=p):
        eps = np.asarray(eps0)
        best = max(best, abs(float(eps @ (K - L) @ eps)) / p**2)
    return best


def d_metric(Kp, Km, Lp, Lm) -> float:
    return 2 * max(q_metric(Kp, Lp), q_metric(Km, Lm))


def trust_value(A: float, B: float) -> tuple[float, float]:
    if A == 0:
        boundary = 0.5 + B / 8
    else:
        boundary = float("inf")

    def objective(t):
        return 0.5 + t / 2 + A / (4 * t) + B / (4 * (t + 2))

    result = minimize_scalar(
        objective, bounds=(1e-12, max(4.0, 2 + sqrt(A + B))), method="bounded"
    )
    if boundary <= result.fun:
        return boundary, 0.0
    return float(result.fun), float(result.x)


def original_alpha_value(mu: float, g: float, h: float, sigma: int) -> float:
    def objective(alpha):
        return alpha + mu**2 * (2 * alpha * g + sigma * h) / (
            2 * (4 * alpha**2 - 1)
        )

    result = minimize_scalar(objective, bounds=(0.50000001, 20.0), method="bounded")
    return float(result.fun)


def random_psd(p: int, rng: np.random.Generator) -> np.ndarray:
    X = rng.normal(size=(p, p))
    K = X @ X.T
    K *= rng.uniform(0.1, 1.0) * p / np.trace(K)
    return K


def verify_sector_metric() -> int:
    rng = np.random.default_rng(260817)
    checks = 0
    for p in range(2, 7):
        Kp, Km, Lp, Lm = (random_psd(p, rng) for _ in range(4))
        G, R = Kp + Km, Kp - Km
        G2, R2 = Lp + Lm, Lp - Lm
        direct = 0.0
        for eps0 in product((-1.0, 1.0), repeat=p):
            eps = np.asarray(eps0)
            for sign in (-1, 1):
                direct = max(
                    direct,
                    abs(float(eps @ ((G - G2) + sign * (R - R2)) @ eps))
                    / p**2,
                )
        assert abs(direct - d_metric(Kp, Km, Lp, Lm)) < 1e-10
        checks += 1
    return checks


def verify_spectral_truncation() -> int:
    rng = np.random.default_rng(170826)
    checks = 0
    for p in (4, 7, 10):
        for eta in (0.2, 0.4, 0.8):
            K = random_psd(p, rng)
            vals, vecs = np.linalg.eigh(K)
            keep = vals > eta * p / 4
            low = (vecs[:, ~keep] * vals[~keep]) @ vecs[:, ~keep].T
            assert int(keep.sum()) <= int(np.floor(4 / eta))
            assert np.linalg.norm(low, 2) <= eta * p / 4 + 1e-10
            assert q_metric(low, np.zeros_like(low)) <= eta / 4 + 1e-10
            checks += 1
    return checks


def verify_trust_coordinates() -> int:
    rng = np.random.default_rng(81726)
    checks = 0
    for _ in range(100):
        g = rng.uniform(0.2, 20.0)
        h = rng.uniform(-g, g)
        mu = rng.uniform(0.05, 1.5)
        sigma = int(rng.choice((-1, 1)))
        a, b = g + sigma * h, g - sigma * h
        transformed, _ = trust_value(mu**2 * a, mu**2 * b)
        original = original_alpha_value(mu, g, h, sigma)
        assert abs(transformed - original) < 2e-5, (
            transformed,
            original,
            mu,
            g,
            h,
        )
        checks += 1
    return checks


def verify_modulus_and_margin() -> int:
    rng = np.random.default_rng(82617)
    checks = 0
    for _ in range(500):
        A, B, A2, B2 = rng.uniform(0, 8, size=4)
        value, t = trust_value(A, B)
        value2, _ = trust_value(A2, B2)
        bound = sqrt(abs(A - A2) / 2) + abs(B - B2) / 8
        assert abs(value - value2) <= bound + 2e-6
        if t > 1e-5:
            derivative_identity = A / t**2 + B / (t + 2) ** 2
            assert abs(derivative_identity - 2) < 5e-5
        tau = rng.uniform(0.05, 2.0)
        criterion = A / tau**2 + B / (tau + 2) ** 2 >= 2
        assert criterion == (t >= tau - 2e-5)
        checks += 1

    for A in (0.01, 0.2, 1.0, 5.0):
        value, t = trust_value(A, 0.0)
        assert abs(value - (0.5 + sqrt(A / 2))) < 2e-6
        assert abs(t - sqrt(A / 2)) < 2e-5
        checks += 1
    return checks


def verify_response_bound() -> int:
    rng = np.random.default_rng(61728)
    checks = 0
    for p in range(2, 7):
        for _ in range(20):
            Kp, Km, Lp, Lm = (random_psd(p, rng) for _ in range(4))
            eta = d_metric(Kp, Km, Lp, Lm)
            mu = rng.uniform(0.05, 1 / p)
            c = mu * p

            def response(Ap, Am):
                best = -float("inf")
                for eps0 in product((-1.0, 1.0), repeat=p):
                    eps = np.asarray(eps0)
                    for sigma in (-1, 1):
                        active = Ap if sigma == 1 else Am
                        passive = Am if sigma == 1 else Ap
                        a = 2 * float(eps @ active @ eps)
                        b = 2 * float(eps @ passive @ eps)
                        value, _ = trust_value(mu**2 * a, mu**2 * b)
                        best = max(best, value)
                return best

            diff = abs(response(Kp, Km) - response(Lp, Lm))
            bound = c * sqrt(eta / 2) + c**2 * eta / 8
            assert diff <= bound + 2e-6
            checks += 1
    return checks


def verify_hard_example() -> int:
    checks = 0
    for p in (2, 5, 9):
        for eta in (0.05, 0.25, 0.8):
            u = np.ones(p) / sqrt(p)
            zero = np.zeros((p, p))
            Kp = eta * p / 2 * np.outer(u, u)
            distance = d_metric(Kp, zero, zero, zero)
            assert abs(distance - eta) < 1e-10
            mu = 1 / p
            A = mu**2 * eta * p**2
            value, _ = trust_value(A, 0.0)
            assert abs(value - 0.5 - sqrt(eta / 2)) < 2e-6
            checks += 1
    return checks


def main() -> None:
    checks = 0
    checks += verify_sector_metric()
    checks += verify_spectral_truncation()
    checks += verify_trust_coordinates()
    checks += verify_modulus_and_margin()
    checks += verify_response_bound()
    checks += verify_hard_example()
    print(f"cross-Gram response metric-entropy checks passed: {checks}")


if __name__ == "__main__":
    main()
