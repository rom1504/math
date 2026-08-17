#!/usr/bin/env python3
"""Exact/numerical checks for the bounded-operator rank barrier."""

from __future__ import annotations

import itertools
import math

import numpy as np


def check_matrix(matrix: np.ndarray) -> None:
    n = matrix.shape[0]
    singular = np.linalg.svd(matrix.astype(float), compute_uv=False)
    c = singular[0] / math.sqrt(n)
    for epsilon in (0.0, 0.2, 0.5, 0.8):
        if epsilon >= 1.0 or c * c <= epsilon * epsilon:
            continue
        visible = int(np.sum(singular > epsilon * math.sqrt(n) + 1e-10))
        lower = n * (1.0 - epsilon**2) / (c**2 - epsilon**2)
        assert visible + 1e-8 >= lower
    rank = int(np.linalg.matrix_rank(matrix.astype(float), tol=1e-9))
    assert singular[0] + 1e-9 >= n / math.sqrt(rank)


def verify() -> None:
    # Exhaust all 4 by 4 sign bridges, then check the two equality models.
    for bits in itertools.product((-1, 1), repeat=16):
        check_matrix(np.array(bits, dtype=int).reshape(4, 4))

    hadamard = np.array(
        [[1, 1, 1, 1], [1, -1, 1, -1],
         [1, 1, -1, -1], [1, -1, -1, 1]],
        dtype=int,
    )
    singular = np.linalg.svd(hadamard.astype(float), compute_uv=False)
    assert np.allclose(singular, 2.0)
    check_matrix(hadamard)
    check_matrix(np.ones((12, 12), dtype=int))
    print("bounded-operator rank barrier checks passed")


if __name__ == "__main__":
    verify()
