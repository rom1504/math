#!/usr/bin/env python3
"""Exact checks for the regular-Hadamard Boolean--spherical trust gap."""

from __future__ import annotations

from itertools import product
from math import sqrt

import numpy as np

from verify_bounded_cap_contextual_metric_compiler import build


V0 = np.asarray(
    (
        -1, -1, -1, 1,
        -1, -1, 1, -1,
        1, -1, 1, 1,
        -1, 1, 1, 1,
    ),
    dtype=np.int64,
)


def boolean_response(H: np.ndarray, a: np.ndarray, b: np.ndarray, m: int) -> int:
    best = 0
    trace = int(np.trace(H))
    for u0 in product((-1, 1), repeat=len(H)):
        u = np.asarray(u0, dtype=np.int64)
        child = (int(u @ H @ u) - trace) // 2
        value = abs(child) + m * (abs(int(a @ u)) + abs(int(b @ u)))
        best = max(best, value)
    return best


def verify_base_exact() -> int:
    r, n, h_list, _ = build(2)
    H = np.asarray(h_list, dtype=np.int64)
    one = np.ones(n, dtype=np.int64)
    assert r == 4 and n == 16
    assert np.trace(H) == 0
    assert np.array_equal(H @ H, r * r * np.eye(n, dtype=np.int64))
    assert np.array_equal(H @ one, r * one)
    assert np.array_equal(H @ V0, r * V0)
    assert int(one @ V0) == 0

    checks = 6
    for m in (r // 2, r):
        boolean = boolean_response(H, one, V0, m)
        expected_boolean = r * n / 2 + m * n
        spherical = r * n / 2 + sqrt(2) * m * n
        assert boolean == expected_boolean
        assert abs((spherical - boolean) - (sqrt(2) - 1) * m * n) < 1e-10
        if m == r // 2:
            assert boolean == 64
            assert abs(spherical - 32 - 32 * sqrt(2)) < 1e-10
        else:
            assert boolean == 96
            assert abs(spherical - 32 - 64 * sqrt(2)) < 1e-10
        checks += 4
    return checks


def verify_tensor_family() -> int:
    r0, n0, h_list, _ = build(2)
    H0 = np.asarray(h_list, dtype=np.int64)
    checks = 0
    H = np.asarray([[1]], dtype=np.int64)
    for j in range(1, 4):
        H = np.kron(H, H0)
        n = n0**j
        r = r0**j
        a = np.ones(n, dtype=np.int64)
        b = np.kron(V0, np.ones(n0 ** (j - 1), dtype=np.int64))
        assert set(a.tolist()) == {1}
        assert set(b.tolist()) <= {-1, 1}
        assert np.array_equal(H @ a, r * a)
        assert np.array_equal(H @ b, r * b)
        assert int(a @ b) == 0
        for eps1, eps2 in product((-1, 1), repeat=2):
            z = eps1 * a + eps2 * b
            assert int(np.sum(np.abs(z))) == n
            assert int(z @ z) == 2 * n
            checks += 2
        m = r // 2
        boolean_formula = r * n / 2 + m * n
        spherical_formula = r * n / 2 + sqrt(2) * m * n
        assert abs(
            (spherical_formula - boolean_formula) / n**1.5
            - (sqrt(2) - 1) / 2
        ) < 1e-12
        total_order = n + 2 * m
        observed_total = (spherical_formula - boolean_formula) / total_order**1.5
        expected_total = (sqrt(2) - 1) / 2 * (n / total_order) ** 1.5
        assert abs(observed_total - expected_total) < 1e-12
        checks += 7
    return checks


def main() -> None:
    checks = verify_base_exact() + verify_tensor_family()
    print(f"regular-Hadamard Boolean--spherical gap checks passed: {checks}")


if __name__ == "__main__":
    main()
