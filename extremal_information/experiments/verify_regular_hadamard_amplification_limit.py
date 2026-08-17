#!/usr/bin/env python3
"""Exact finite checks for regular-Hadamard extremal amplification."""

from __future__ import annotations

from itertools import product

import numpy as np


def quadratic_max(matrix: np.ndarray, absolute: bool = False) -> int:
    best = -10**18
    for values in product((-1, 1), repeat=len(matrix)):
        x = np.asarray(values, dtype=np.int64)
        value = int(x @ matrix @ x)
        if absolute:
            value = abs(value)
        best = max(best, value)
    return best


def verify() -> None:
    h = np.asarray(
        [[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]],
        dtype=np.int64,
    )
    u = np.asarray([1, 1, 1, -1], dtype=np.int64)
    assert np.array_equal(h @ h, 4 * np.eye(4, dtype=np.int64))
    assert np.array_equal(h @ u, 2 * u)
    assert int(u @ h @ u) == 8

    matrices = (
        np.asarray([[1, 1], [1, -1]], dtype=np.int64),
        np.asarray([[2, -1], [-1, 0]], dtype=np.int64),
        np.asarray([[1, -2, 1], [-2, 0, 1], [1, 1, -1]], dtype=np.int64),
    )
    checks = 3
    for b in matrices:
        d = len(b)
        b1 = np.kron(b, h)
        for absolute in (False, True):
            q0 = quadratic_max(b, absolute) / (2 * d ** 1.5)
            q1 = quadratic_max(b1, absolute) / (2 * (4 * d) ** 1.5)
            assert q1 + 1e-12 >= q0
            assert q1 <= np.linalg.norm(b, 2) / (2 * np.sqrt(d)) + 1e-12
            checks += 2

        # Every level-one Boolean witness embeds at level two with exactly
        # the same normalized pair correlations and quadratic value.
        rng = np.random.default_rng(20260817 + d)
        for _ in range(16):
            x = rng.choice(np.asarray([-1, 1], dtype=np.int64), size=4 * d)
            xu = np.kron(x, u)
            left = int(x @ b1 @ x)
            right = int(xu @ np.kron(b1, h) @ xu)
            assert right == 8 * left
            checks += 1

    # A full sign template with zero trace amplifies to a valid hollow sign
    # matrix without changing its Boolean quadratic energy.
    b = np.asarray(
        [[1, 1, -1, 1], [1, 1, 1, -1], [-1, 1, -1, 1], [1, -1, 1, -1]],
        dtype=np.int64,
    )
    assert int(np.trace(b)) == 0
    amplified = np.kron(b, h)
    hollow = amplified.copy()
    np.fill_diagonal(hollow, 0)
    assert set(int(v) for v in hollow[~np.eye(len(hollow), dtype=bool)]) == {-1, 1}
    rng = np.random.default_rng(42)
    for _ in range(64):
        x = rng.choice(np.asarray([-1, 1], dtype=np.int64), size=len(amplified))
        assert int(x @ amplified @ x) == int(x @ hollow @ x)
        checks += 1

    print(f"regular-Hadamard amplification checks passed: {checks}")


if __name__ == "__main__":
    verify()
