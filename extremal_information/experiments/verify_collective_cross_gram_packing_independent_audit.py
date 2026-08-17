#!/usr/bin/env python3
"""Independent checks for collective cross-Gram packing and response modulus.

This supplements the canonical verifier with exact projective-ball counts,
pointwise SA.19/CP.11 normalization, and an explicit common-involution
realization of the hard pair.
"""

from __future__ import annotations

import itertools
from math import comb, sqrt

import numpy as np

from verify_collective_cross_gram_packing import collective_distance, psi


def sign_words(p: int) -> list[np.ndarray]:
    return [np.asarray(word, dtype=int) for word in itertools.product((-1, 1), repeat=p)]


def projective_distance(left: np.ndarray, right: np.ndarray) -> int:
    distance = int(np.count_nonzero(left != right))
    return min(distance, len(left) - distance)


def check_metric_factor_and_projective_balls() -> int:
    checks = 0
    for p in range(2, 7):
        words = [word for word in sign_words(p) if word[0] == 1]
        for s in words:
            for t in words:
                h = int(np.count_nonzero(s != t))
                g = np.outer(s, s)
                gp = np.outer(t, t)
                observed = collective_distance(g, g, gp, gp)
                expected = 8 * h * (p - h) / p**2
                assert abs(observed - expected) < 1e-12
                assert abs(expected - 2 * (1 - (float(s @ t) / p) ** 2)) < 1e-12
                checks += 1

        center = np.ones(p, dtype=int)
        for radius in range(1, p // 2 + 1):
            observed_ball = sum(projective_distance(center, word) < radius for word in words)
            expected_ball = sum(comb(p, level) for level in range(radius))
            assert observed_ball == expected_ball
            checks += 1
    return checks


def check_sa19_normalization() -> int:
    checks = 0
    rng = np.random.default_rng(170826)
    for p in (2, 5, 11):
        for sigma in (-1, 1):
            for _ in range(100):
                a, b = rng.random(2)
                kappa = float(rng.uniform(0.05, 3.0))
                t = float(rng.uniform(1e-4, 5.0))
                alpha = (1 + t) / 2
                dangerous = p * p * a
                safe = p * p * b
                if sigma == 1:
                    k_plus, k_minus = dangerous, safe
                else:
                    k_plus, k_minus = safe, dangerous
                g = k_plus + k_minus
                h = k_plus - k_minus
                sa19 = alpha + (kappa / p) ** 2 * (2 * alpha * g + sigma * h) / (
                    2 * (4 * alpha * alpha - 1)
                )
                cp11 = (1 + t) / 2 + kappa**2 / 2 * (a / t + b / (t + 2))
                assert abs(sa19 - cp11) < 2e-11
                checks += 1
    return checks


def sectors(g: np.ndarray, r: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    return (g + r) / 2, (g - r) / 2


def gram_rayleigh(ports: np.ndarray, involution: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    n = ports.shape[1]
    return ports @ ports.T / n, ports @ involution @ ports.T / n


def check_hard_realization() -> int:
    p, n = 2, 4
    s = np.asarray((1, 1), dtype=float)
    t = np.asarray((1, -1), dtype=float)
    involution = np.diag((1.0, 1.0, -1.0, -1.0))
    checks = 0
    for c in (1e-6, 0.01, 0.25, 0.75):
        ports_a = np.zeros((p, n))
        ports_b = np.zeros((p, n))
        for i in range(p):
            ports_a[i, 0] = sqrt(n * c) * s[i]
            ports_b[i, 0] = sqrt(n * c) * t[i]
            ports_a[i, 2 + i] = sqrt(n * (1 - c))
            ports_b[i, 2 + i] = sqrt(n * (1 - c))

        ga, ra = gram_rayleigh(ports_a, involution)
        gb, rb = gram_rayleigh(ports_b, involution)
        kap, kam = sectors(ga, ra)
        kbp, kbm = sectors(gb, rb)
        assert np.allclose(kap, c * np.outer(s, s))
        assert np.allclose(kbp, c * np.outer(t, t))
        assert np.allclose(kam, (1 - c) * np.eye(p))
        assert np.allclose(kbm, kam)
        assert np.allclose(np.diag(ga), 1)
        assert np.allclose(np.diag(gb), 1)
        assert np.allclose(np.diag(ra), 2 * c - 1)
        assert np.allclose(np.diag(rb), 2 * c - 1)
        assert abs(collective_distance(ga, ra, gb, rb) - 2 * c) < 1e-12

        a0 = float(t @ kap @ t) / p**2
        a1 = float(t @ kbp @ t) / p**2
        b0 = float(t @ kam @ t) / p**2
        assert abs(a0) < 1e-12 and abs(a1 - c) < 1e-12
        assert abs(b0 - (1 - c) / 2) < 1e-12
        checks += 1
    return checks


def check_unrestricted_kappa_sharp_variant() -> int:
    """A rank-one safe sector gives b=0 and an exact sqrt gap for any kappa."""

    s = np.asarray((1, 1), dtype=float)
    t = np.asarray((1, -1), dtype=float)
    p = 2
    checks = 0
    for kappa in (0.2, 1.0, 4.0):
        for c in (1e-8, 1e-4, 0.1):
            k_plus_a = c * np.outer(s, s)
            k_plus_b = c * np.outer(t, t)
            k_minus = (1 - c) * np.outer(s, s)
            a0 = float(t @ k_plus_a @ t) / p**2
            a1 = float(t @ k_plus_b @ t) / p**2
            b = float(t @ k_minus @ t) / p**2
            assert abs(a0) < 1e-12 and abs(a1 - c) < 1e-12 and abs(b) < 1e-12
            gap = psi(kappa, a1, b) - psi(kappa, a0, b)
            assert abs(gap - kappa * sqrt(c)) < 2e-11
            checks += 1
    return checks


def main() -> None:
    metric_checks = check_metric_factor_and_projective_balls()
    normalization_checks = check_sa19_normalization()
    hard_checks = check_hard_realization()
    variant_checks = check_unrestricted_kappa_sharp_variant()
    print(
        "collective cross-Gram independent audit: PASS",
        f"metric_ball={metric_checks}",
        f"normalization={normalization_checks}",
        f"hard_realization={hard_checks}",
        f"sharp_variant={variant_checks}",
    )


if __name__ == "__main__":
    main()
