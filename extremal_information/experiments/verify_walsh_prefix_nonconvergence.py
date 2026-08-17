#!/usr/bin/env python3
"""Exact certificate for two separated prefix phases of infinite Walsh."""

from __future__ import annotations

import math

import numpy as np


H = np.asarray(
    [[1, 1, 1, 1], [1, -1, 1, -1],
     [1, 1, -1, -1], [1, -1, -1, 1]],
    dtype=np.int64,
)
U = np.asarray([1, 1, 1, -1], dtype=np.int64)
Z = np.asarray(
    [
        1, 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, -1,
        1, -1, -1, -1, 1, 1, 1, -1, 1, -1, 1, 1, 1, 1, -1, -1,
        1, 1, 1, 1, -1, 1, -1, 1, -1, 1, 1, -1, -1, -1, 1, -1,
    ],
    dtype=np.int64,
)


def verify() -> None:
    checks = 0
    assert np.array_equal(H @ H, 4 * np.eye(4, dtype=np.int64))
    assert np.array_equal(H @ U, 2 * U)
    assert int(np.trace(H)) == 0
    assert set(map(int, Z)) == {-1, 1}
    checks += 4

    powers = [np.asarray([[1]], dtype=np.int64)]
    for _ in range(4):
        powers.append(np.kron(H, powers[-1]))
    for r in range(1, 4):
        assert np.array_equal(powers[r + 1][: 4**r, : 4**r], powers[r])
        pole = U.copy()
        for _ in range(r - 1):
            pole = np.kron(U, pole)
        assert np.array_equal(powers[r] @ pole, (2**r) * pole)
        assert int(pole @ powers[r] @ pole) == (4**r) ** 1.5
        assert int(np.trace(powers[r])) == 0
        checks += 4

    # Every fixed base-four prefix phase factors as one fixed outer prefix
    # tensor a deeper Walsh power.  This is the algebraic input to WP.2.
    prefix_checks = 0
    for r in range(1, 4):
        ambient = powers[r + 1]
        for k in range(r + 1):
            outer = powers[k + 1]
            for p in range(4**k, 4 ** (k + 1) + 1):
                size = p * 4 ** (r - k)
                expected_prefix = np.kron(
                    outer[:p, :p], powers[r - k]
                )
                assert np.array_equal(
                    ambient[:size, :size], expected_prefix
                )
                prefix_checks += 1
    checks += prefix_checks

    B = H[:3, :3]
    base = np.kron(B, powers[2])
    assert len(base) == 48 == len(Z)
    assert int(Z @ base @ Z) == 356
    assert int(np.trace(base)) == 0
    checks += 3

    for r in range(1, 4):
        assert np.array_equal(
            powers[r + 1][: 3 * 4**r, : 3 * 4**r],
            np.kron(B, powers[r]),
        )
        checks += 1

    witness = Z.copy()
    matrix = base.copy()
    expected = 356
    for _ in range(2):
        witness = np.kron(witness, U)
        matrix = np.kron(matrix, H)
        expected *= 8
        assert int(witness @ matrix @ witness) == expected
        checks += 1

    ratio = 89 / (96 * math.sqrt(3))
    assert ratio > 0.5
    assert 89 * 89 > 48 * 48 * 3
    checks += 2

    print(
        "Walsh prefix nonconvergence checks passed: "
        f"{checks}; high={ratio:.15f}; gap={ratio - 0.5:.15f}"
    )


if __name__ == "__main__":
    verify()
