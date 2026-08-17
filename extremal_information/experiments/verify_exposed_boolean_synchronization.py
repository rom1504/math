#!/usr/bin/env python3
"""Checks for exposed Boolean synchronization and the Walsh close-pole family."""

from __future__ import annotations

from itertools import product
from math import gamma, pi, sqrt

import numpy as np


def bits(d: int):
    return list(product((0, 1), repeat=d))


def walsh_matrix(dim: int) -> np.ndarray:
    vv = bits(dim)
    return np.asarray(
        [[(-1) ** (sum(a * b for a, b in zip(x, y)) % 2) for y in vv] for x in vv],
        dtype=np.int64,
    )


def pairing(x):
    out = list(x)
    for i in range(0, len(x), 2):
        out[i], out[i + 1] = x[i + 1], x[i]
    return tuple(out)


def construction(d: int):
    q = 2**d
    small = walsh_matrix(d)
    W = np.kron(small, small)
    vv = bits(d)
    y0 = np.asarray(
        [(-1) ** (sum(a * b for a, b in zip(x, z)) % 2) for x in vv for z in vv],
        dtype=np.int64,
    )
    L = np.asarray(
        [1 if z == pairing(x) else 0 for x in vv for z in vv], dtype=np.int64
    )
    y1 = y0 - 2 * L
    H = y0[:, None] * W * y0[None, :]
    a = np.ones(q * q, dtype=np.int64)
    b = y0 * y1
    return q, W, H, y0, y1, L, a, b


def boolean_response(H, a, b, width):
    best = 0
    trace = int(np.trace(H))
    for x0 in product((-1, 1), repeat=len(H)):
        x = np.asarray(x0, dtype=np.int64)
        child = abs((int(x @ H @ x) - trace) // 2)
        field = width * (abs(int(a @ x)) + abs(int(b @ x)))
        best = max(best, child + field)
    return best


def verify_walsh_family() -> int:
    checks = 0
    for d in (2, 4):
        q, W, H, y0, y1, L, a, b = construction(d)
        n = q * q
        assert np.array_equal(W @ W, n * np.eye(n, dtype=np.int64))
        assert np.array_equal(W @ y0, q * y0)
        assert np.array_equal(W @ L, q * L)
        assert np.array_equal(W @ y1, q * y1)
        assert int(np.sum(L)) == q
        assert np.all(y0[L.astype(bool)] == 1)
        assert set(y1.tolist()) <= {-1, 1}
        assert np.array_equal(H @ H, n * np.eye(n, dtype=np.int64))
        assert np.trace(H) == 0
        assert np.array_equal(H @ a, q * a)
        assert np.array_equal(H @ b, q * b)
        assert int(a @ b) == n - 2 * q

        rho = 1 - 2 / q
        u = (a + b) / sqrt(2 * (1 + rho))
        assert abs(float(u @ u) - n) < 1e-9
        flatness = float(np.sum(np.abs(u))) / n
        assert abs(flatness - sqrt(1 - 1 / q)) < 1e-12

        width = q // 2
        boolean_formula = q * n / 2 + width * n * (1 + rho)
        spherical_formula = q * n / 2 + width * n * sqrt(2 * (1 + rho))
        normalized_gap = (spherical_formula - boolean_formula) / (q * n)
        expected_gap = sqrt(1 - 1 / q) - (1 - 1 / q)
        assert abs(normalized_gap - expected_gap) < 1e-12
        if d == 2:
            assert boolean_response(H, a, b, width) == boolean_formula
        checks += 17
    return checks


def verify_rounding_bound() -> int:
    rng = np.random.default_rng(817260)
    checks = 0
    for n in (8, 16):
        # Symmetric involution with operator norm r=1.
        Q, _ = np.linalg.qr(rng.normal(size=(n, n)))
        signs = np.diag(rng.choice((-1.0, 1.0), size=n))
        H = Q @ signs @ Q.T
        ports = rng.choice((-1.0, 1.0), size=(3, n))
        m = 0.4
        c = m * len(ports)
        for _ in range(100):
            u = rng.normal(size=n)
            u *= sqrt(n) / np.linalg.norm(u)
            sigma = int(rng.choice((-1, 1)))
            eps = rng.choice((-1.0, 1.0), size=len(ports))
            z = eps @ ports
            x = np.where(u >= 0, 1.0, -1.0)
            sphere = sigma * float(u @ H @ u) / 2 + m * float(z @ u)
            cube = sigma * float(x @ H @ x) / 2 + m * float(z @ x)
            phi = 1 - np.sum(np.abs(u)) / n
            bound = (1 + c) * n * sqrt(2 * phi)
            assert sphere - cube <= bound + 1e-9
            checks += 1
    return checks


def verify_subspace_constant() -> int:
    checks = 0
    for d in range(1, 9):
        gamma_d = sqrt(d) * gamma(d / 2) / (sqrt(pi) * gamma((d + 1) / 2))
        if d == 1:
            assert abs(gamma_d - 1) < 1e-12
        else:
            assert gamma_d < 1
        checks += 1
    return checks


def main() -> None:
    checks = verify_walsh_family()
    checks += verify_rounding_bound()
    checks += verify_subspace_constant()
    print(f"exposed Boolean synchronization checks passed: {checks}")


if __name__ == "__main__":
    main()
