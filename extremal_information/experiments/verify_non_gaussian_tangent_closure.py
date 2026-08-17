#!/usr/bin/env python3
"""Numerical checks for power-roof and tangent-shape closure."""

from __future__ import annotations

import math

import numpy as np


def effective(a: float, b: float, p: float) -> float:
    return (a ** (-1.0 / (p - 1.0)) + b ** (-1.0 / (p - 1.0))) ** (-(p - 1.0))


def verify_roof() -> None:
    grid = np.linspace(-2.0, 3.0, 500001)
    for p in (1.5, 2.0, 4.0):
        for a, b, z in ((0.7, 1.9, 1.3), (2.0, 0.4, -0.8)):
            observed = float(np.min(a * np.abs(grid) ** p + b * np.abs(z - grid) ** p))
            target = effective(a, b, p) * abs(z) ** p
            assert abs(observed - target) < 2e-8


def standardized_kurtosis(p: float) -> float:
    # E X^(2j)=a^(-2j/p) Gamma((2j+1)/p)/Gamma(1/p).
    return math.gamma(5.0 / p) * math.gamma(1.0 / p) / math.gamma(3.0 / p) ** 2


def verify_tangent_obstruction() -> None:
    # The self-convolution halves excess kurtosis.  It can remain in the same
    # scale family only when excess kurtosis is zero; among these test powers
    # that occurs exactly at p=2.
    assert abs(standardized_kurtosis(2.0) - 3.0) < 1e-12
    for p in (1.25, 1.5, 3.0, 4.0, 6.0):
        excess = standardized_kurtosis(p) - 3.0
        assert abs(excess) > 1e-3
        assert abs(excess / 2.0 - excess) > 1e-3

    # Central lattice mass has exponent 1-1/p.
    p = 4.0
    a = 0.8
    b = 1.1
    integral = 2.0 * math.gamma(1.0 / p) / (p * (a + b) ** (1.0 / p))
    errors = []
    for n in (100, 300, 1000, 3000):
        radius = int(12 * n ** (1.0 - 1.0 / p))
        k = np.arange(-radius, radius + 1, dtype=float)
        mass = float(np.sum(np.exp(-(a + b) * np.abs(k) ** p / n ** (p - 1.0))))
        scaled = mass / n ** (1.0 - 1.0 / p)
        errors.append(abs(scaled - integral))
    # For this analytic even profile Poisson-summation errors are already at
    # floating precision at these sizes; monotonicity is not numerically
    # meaningful.
    assert max(errors) < 2e-4


def verify() -> None:
    verify_roof()
    verify_tangent_obstruction()
    print("non-Gaussian tangent closure checks passed")


if __name__ == "__main__":
    verify()
