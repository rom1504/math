#!/usr/bin/env python3
"""Finite wind tunnel for the Morse tangent-mass composition law."""

from __future__ import annotations

import math

import numpy as np


def log_binomial(n: int, k: int) -> float:
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def logsumexp(values: list[float]) -> float:
    top = max(values)
    return top + math.log(sum(math.exp(value - top) for value in values))


def verify() -> None:
    checks = 0
    # Vandermonde: pointwise maximum misses a sqrt(n) saddle mass.
    gaps = []
    for n in (20, 40, 80, 160, 320, 640):
        target = n  # p=1/2 in a total size 2n convolution
        terms = [
            log_binomial(n, k) + log_binomial(n, target - k)
            for k in range(n + 1)
        ]
        total = logsumexp(terms)
        exact = log_binomial(2 * n, target)
        assert abs(total - exact) < 2e-10
        gap = total - max(terms)
        gaps.append(gap / math.log(n))
        amplitude_ratio = math.exp(
            exact - (2 * n * math.log(2) - 0.5 * math.log(math.pi * n))
        )
        assert abs(amplitude_ratio - 1.0) < 0.02
        checks += 1
    assert gaps[-1] > 0.45
    assert abs(gaps[-1] - 0.5) < abs(gaps[0] - 0.5)
    checks += 2

    # A quartic saddle has n^(3/4), falsifying a universal d/2 correction.
    quartic_slopes = []
    for n in (100, 400, 1600, 6400, 25600):
        logs = [-(k**4) / (n**3) for k in range(-n, n + 1)]
        quartic_slopes.append(logsumexp(logs) / math.log(n))
        checks += 1
    assert quartic_slopes[-1] > 0.70
    assert abs(quartic_slopes[-1] - 0.75) < abs(quartic_slopes[0] - 0.75)
    quartic_constant = math.gamma(0.25) / 2
    quartic_sum = math.exp(logsumexp(
        [-(k**4) / (25600**3) for k in range(-25600, 25601)]
    ))
    assert abs(quartic_sum / 25600**0.75 - quartic_constant) < 0.04
    checks += 2

    # A genuinely two-dimensional, off-centre Gaussian tangent convolution.
    P = np.asarray([[2.0, 0.35], [0.35, 1.4]])
    Q = np.asarray([[1.1, -0.2], [-0.2, 1.8]])
    mu = np.asarray([0.2, -0.1])
    nu = np.asarray([-0.15, 0.25])
    z0 = np.asarray([0.45, -0.2])
    a, b = 1.3, 0.8
    R = np.linalg.inv(np.linalg.inv(P) + np.linalg.inv(Q))
    a_out = (2 * math.pi) * a * b / math.sqrt(np.linalg.det(P + Q))
    gaussian_ratios = []
    for n in (25, 64, 121):
        ell = np.rint(n * z0).astype(int)
        z = ell / n
        saddle = np.linalg.solve(P + Q, P @ mu + Q @ (z - nu))
        radius = int(math.ceil(8 * math.sqrt(n)))
        centre = np.rint(n * saddle).astype(int)
        g0 = np.arange(centre[0] - radius, centre[0] + radius + 1)
        g1 = np.arange(centre[1] - radius, centre[1] + radius + 1)
        K0, K1 = np.meshgrid(g0, g1, indexing="ij")
        points = np.stack([K0.ravel(), K1.ravel()], axis=1)
        left = points - n * mu
        right = ell - points - n * nu
        exponents = -(
            np.einsum("bi,ij,bj->b", left, P, left)
            + np.einsum("bi,ij,bj->b", right, Q, right)
        ) / (2 * n)
        top = float(np.max(exponents))
        actual = a * b * math.exp(top) * float(
            np.exp(exponents - top).sum()
        )
        leading = -0.5 * n * (z - mu - nu) @ R @ (z - mu - nu)
        predicted = n * a_out * math.exp(leading)
        gaussian_ratios.append(actual / predicted)
        checks += 1
    assert max(abs(ratio - 1) for ratio in gaussian_ratios) < 2e-8
    checks += 1

    # Restricting a one-dimensional Gaussian to 2Z divides tangent density
    # by the covolume two.
    n = 400
    even_sum = sum(
        math.exp(-(k * k) / (2 * n))
        for k in range(-8 * int(math.sqrt(n)), 8 * int(math.sqrt(n)) + 1, 2)
    )
    assert abs(even_sum / (math.sqrt(2 * math.pi * n) / 2) - 1) < 1e-10
    checks += 1

    print(
        "Morse tangent-mass checks passed: "
        f"{checks}; Vandermonde slopes={gaps}; quartic={quartic_slopes}; "
        f"Gaussian ratios={gaussian_ratios}"
    )


if __name__ == "__main__":
    verify()
