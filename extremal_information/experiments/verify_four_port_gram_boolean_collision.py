#!/usr/bin/env python3
"""Exact checks for the four-port equal-Gram Boolean support collision."""

from __future__ import annotations

from itertools import product

import numpy as np


def systems(L: int) -> tuple[np.ndarray, np.ndarray]:
    rows = np.asarray(list(product((-1, 1), repeat=4)), dtype=np.int64)
    even = rows[np.prod(rows, axis=1) == 1]
    U = np.repeat(rows, L, axis=0)
    V = np.repeat(even, 2 * L, axis=0)
    return U, V


def endpoint_response(W: np.ndarray) -> int:
    eps = np.asarray(list(product((-1, 1), repeat=4)), dtype=np.int64)
    return int(np.max(np.sum(np.abs(eps @ W.T), axis=1)))


def old_spin_response(W: np.ndarray) -> int:
    xs = np.asarray(list(product((-1, 1), repeat=len(W))), dtype=np.int64)
    return int(np.max(np.sum(np.abs(xs @ W), axis=1)))


def main() -> None:
    checks = 0
    for L in (1, 2, 5, 9):
        U, V = systems(L)
        n = 16 * L
        assert np.array_equal(U.T @ U, n * np.eye(4, dtype=np.int64))
        assert np.array_equal(V.T @ V, n * np.eye(4, dtype=np.int64))
        assert endpoint_response(U) == 3 * n // 2
        assert endpoint_response(V) == 2 * n
        if L == 1:
            assert old_spin_response(U) == endpoint_response(U)
            assert old_spin_response(V) == endpoint_response(V)
        checks += 6

    U, V = systems(1)
    eps = np.asarray(list(product((-1, 1), repeat=4)), dtype=np.int64)
    even_values = sorted(set(np.sum(np.abs(eps @ V.T), axis=1).tolist()))
    assert even_values == [16, 32]
    assert np.allclose(np.linalg.svd(U, compute_uv=False), 4)
    assert np.allclose(np.linalg.svd(V, compute_uv=False), 4)
    checks += 3
    print(f"four-port Gram collision checks passed: {checks}")


if __name__ == "__main__":
    main()
