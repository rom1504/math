#!/usr/bin/env python3
"""Exact finite checks for the planted random-partition identities."""

from __future__ import annotations

import itertools
import math
import random

import numpy as np


def trial(n: int, k: int, seed: int) -> None:
    rng = random.Random(seed)
    upper = np.array(
        [[0 if i == j else rng.choice((-1, 1)) for j in range(n)]
         for i in range(n)],
        dtype=float,
    )
    A = np.triu(upper, 1)
    A = A + A.T
    x = np.array([rng.choice((-1, 1)) for _ in range(n)], dtype=float)
    g = np.array([rng.randrange(k) for _ in range(n)], dtype=int)
    P = np.zeros((n, k))
    P[np.arange(n), g] = 1
    D = np.diag(x)

    ebar = np.ones(k) / k
    Z = np.zeros((n, k))
    left = np.zeros((n, n))
    right = np.zeros((k, k))
    for i in range(n):
        vi = x[i] * A[:, i]
        # Enumerate the expectation over the block of this one summand.
        for a in range(k):
            w = np.eye(k)[a] - ebar
            Zi = np.outer(vi, w)
            left += (Zi @ Zi.T) / k
            right += (Zi.T @ Zi) / k
        Z += np.outer(vi, np.eye(k)[g[i]] - ebar)

    M = A @ D @ P
    assert np.allclose(M, np.outer(A @ x, ebar) + Z)
    assert np.allclose(left, (1 - 1 / k) * (A @ A))
    expected_right = n * (n - 1) * (
        np.eye(k) / k - np.ones((k, k)) / (k * k)
    )
    assert np.allclose(right, expected_right)

    exact_max = max(
        np.linalg.norm(M @ np.array(z, dtype=float)) ** 2
        for z in itertools.product((-1, 1), repeat=k)
    )
    spectral_bound = k * np.linalg.norm(M, 2) ** 2
    assert exact_max <= spectral_bound + 1e-8


def main() -> None:
    for n in range(5, 10):
        for k in range(2, min(5, n) + 1):
            for seed in range(8):
                trial(n, k, 1000 * n + 100 * k + seed)
    print("PASS: planted-partition decomposition, variances, and spectral bound")


if __name__ == "__main__":
    main()
