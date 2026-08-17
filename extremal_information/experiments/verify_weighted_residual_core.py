#!/usr/bin/env python3
"""Exact checks for the bounded-delay residual/core separation."""

from __future__ import annotations

from itertools import product

import numpy as np


TA = np.asarray([[0, -1], [-2, -3]], dtype=np.int64)
TB = np.asarray([[-2, -3], [1, 0]], dtype=np.int64)
P = np.asarray([0, -1], dtype=np.int64)


def mp_row(row: np.ndarray, matrix: np.ndarray) -> np.ndarray:
    return np.max(row[:, None] + matrix, axis=0)


def mp_product(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return np.max(left[:, :, None] + right[None, :, :], axis=1)


def max_cycle_mean(matrix: np.ndarray) -> float:
    # Two-state exact enumeration: the maximizing closed walk is a loop or
    # the two-cycle.  This is enough for the scoped counterexample.
    return max(
        float(matrix[0, 0]),
        float(matrix[1, 1]),
        float(matrix[0, 1] + matrix[1, 0]) / 2,
    )


def threshold_image(matrix: np.ndarray, subset: frozenset[int]) -> frozenset[int]:
    return frozenset(
        j for j in range(2)
        if any(matrix[i, j] >= 0 for i in subset)
    )


def verify() -> None:
    checks = 0
    for matrix in (TA, TB):
        assert np.array_equal(mp_row(P, matrix), P)
        checks += 1

    # Every word through length eight retains p as a zero-eigenprofile and
    # has spectral radius zero.
    for length in range(1, 9):
        for word in product((TA, TB), repeat=length):
            matrix = word[0]
            for next_matrix in word[1:]:
                matrix = mp_product(matrix, next_matrix)
            assert np.array_equal(mp_row(P, matrix), P)
            assert max_cycle_mean(matrix) == 0
            checks += 2

    assert threshold_image(TA, frozenset({0, 1})) == frozenset({0})
    assert threshold_image(TB, frozenset({0, 1})) == frozenset({0, 1})
    checks += 2

    # Descending incoming-image iteration for the single residual context.
    core = frozenset({0, 1})
    while True:
        new_core = (
            core
            & threshold_image(TA, core)
            & threshold_image(TB, core)
        )
        if new_core == core:
            break
        core = new_core
    assert core == frozenset()
    checks += 1

    # Coordinate pins attain the endpoint differences underlying
    # half-oscillation projective distance.
    p = np.asarray([0.0, -1.0, -3.0])
    q = np.asarray([-2.0, 0.0, -1.0])
    diff = p - q
    values = []
    for coordinate in range(3):
        z = np.full(3, -1000.0)
        z[coordinate] = 0.0
        values.append(float(np.max(p + z) - np.max(q + z)))
    assert np.allclose(values, diff)
    assert (max(values) - min(values)) / 2 == 2.0
    checks += 2

    print(f"weighted residual/core checks passed: {checks}")


if __name__ == "__main__":
    verify()
