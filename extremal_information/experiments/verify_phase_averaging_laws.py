#!/usr/bin/env python3
"""Numerical wind tunnel for the automatic phase averaging formulas."""

from __future__ import annotations

import math

import numpy as np


H = 4


def phase(t: np.ndarray | float) -> np.ndarray | float:
    return 0.5 + 0.04 * np.sin(math.pi * np.log(t) / math.log(H)) ** 2


def trap_integral(values: np.ndarray, points: np.ndarray) -> float:
    return float(
        np.sum((values[:-1] + values[1:]) * (points[1:] - points[:-1])) / 2.0
    )


def verify() -> None:
    grid = np.linspace(1.0, float(H), 200001)
    log_target = trap_integral(phase(grid) / grid, grid) / math.log(H)
    full_integral = trap_integral(phase(grid), grid)

    log_errors = []
    cesaro_errors = []
    power_errors = {0.3: [], 1.7: []}
    s = 2.7
    partial_grid = np.linspace(1.0, s, 100001)
    cesaro_target = (
        full_integral / (H - 1)
        + trap_integral(phase(partial_grid), partial_grid)
    ) / s

    for r in (4, 5, 6, 7, 8):
        n_max = int(s * H**r)
        n = np.arange(1, n_max + 1, dtype=np.int64)
        powers = np.floor(np.log(n) / math.log(H)).astype(int)
        # Correct occasional floating logarithm errors at exact powers.
        base = H**powers
        too_large = base > n
        powers[too_large] -= 1
        base = H**powers
        while np.any(H * base <= n):
            mask = H * base <= n
            powers[mask] += 1
            base[mask] *= H
        t = n / base
        # One fixed sequence satisfying the theorem's block-uniform error,
        # rather than a triangular array depending on the cutoff r.
        q = phase(t) + 0.02 / np.sqrt(base)
        log_mean = float(np.sum(q / n) / math.log(n_max))
        ordinary = float(np.mean(q))
        log_errors.append(abs(log_mean - log_target))
        cesaro_errors.append(abs(ordinary - cesaro_target))
        for alpha in power_errors:
            weights = n.astype(float) ** (alpha - 1.0)
            full_weighted = trap_integral(
                grid ** (alpha - 1.0) * phase(grid), grid
            )
            partial_weighted = trap_integral(
                partial_grid ** (alpha - 1.0) * phase(partial_grid),
                partial_grid,
            )
            target = alpha * (
                full_weighted / (H**alpha - 1.0) + partial_weighted
            ) / s**alpha
            observed = float(np.sum(weights * q) / np.sum(weights))
            power_errors[alpha].append(abs(observed - target))

    assert log_errors[-1] < log_errors[0]
    assert cesaro_errors[-1] < cesaro_errors[0]
    # Harmonic-block convergence is only O(1/log N) because early blocks
    # retain a fixed numerator contribution.
    assert log_errors[-1] < 0.06
    assert cesaro_errors[-1] < 0.002
    for errors in power_errors.values():
        assert errors[-1] < errors[0]
        assert errors[-1] < 0.015
    print(
        "phase averaging checks passed; "
        f"log errors={log_errors}; Cesaro errors={cesaro_errors}; "
        f"power errors={power_errors}"
    )


if __name__ == "__main__":
    verify()
