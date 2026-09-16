#!/usr/bin/env python3
"""Exact rational checks for heterogeneous product-noise identities."""

from __future__ import annotations

from fractions import Fraction
from itertools import product


def matvec(a: list[list[int]], x: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sum(a[i][j] * x[j] for j in range(len(x))) for i in range(len(x)))


def quad(a: list[list[int]], x: tuple[int, ...]) -> int:
    ax = matvec(a, x)
    return sum(x[i] * ax[i] for i in range(len(x)))


def deterministic_signing(n: int) -> list[list[int]]:
    a = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            value = 1 if ((i * j + 2 * i + j * j + 3 * j) % 7) < 4 else -1
            a[i][j] = a[j][i] = value
    return a


def complement_matrix(a: list[list[int]], selector: frozenset[int]) -> list[list[int]]:
    n = len(a)
    return [
        [
            0 if i == j else a[i][j] * (1 if i in selector and j in selector else -1)
            for j in range(n)
        ]
        for i in range(n)
    ]


def exact_check(n: int, means: tuple[Fraction, ...]) -> None:
    assert len(means) == n
    a = deterministic_signing(n)
    b = complement_matrix(a, frozenset(range(0, n, 2)))
    base = tuple(1 if (3 * i + i * i) % 5 < 3 else -1 for i in range(n))
    sigma = -1
    c = [
        [sigma * base[i] * b[i][j] * base[j] for j in range(n)]
        for i in range(n)
    ]

    c_rho = tuple(
        sum(Fraction(c[i][j]) * means[j] for j in range(n))
        for i in range(n)
    )
    expected_energy = sum(means[i] * c_rho[i] for i in range(n))
    variances = tuple(1 - r * r for r in means)
    expected_variance = (
        4 * sum(variances[i] * c_rho[i] * c_rho[i] for i in range(n))
        + 4 * sum(
            variances[i] * variances[j]
            for i in range(n)
            for j in range(i + 1, n)
        )
    )

    # A D_base rho and the independent-column variance in (R38.A3).
    mean_spin = tuple(Fraction(base[i]) * means[i] for i in range(n))
    a_mean = tuple(
        sum(Fraction(a[i][j]) * mean_spin[j] for j in range(n))
        for i in range(n)
    )
    expected_row = sum(v * v for v in a_mean) + (n - 1) * sum(variances)

    total = Fraction(0)
    mean_energy = Fraction(0)
    mean_energy2 = Fraction(0)
    mean_row = Fraction(0)
    for z in product((-1, 1), repeat=n):
        probability = Fraction(1)
        for i, value in enumerate(z):
            probability *= (1 + value * means[i]) / 2
        noisy = tuple(base[i] * z[i] for i in range(n))
        energy = sigma * quad(b, noisy)
        row_vector = matvec(a, noisy)
        row = sum(v * v for v in row_vector)
        total += probability
        mean_energy += probability * energy
        mean_energy2 += probability * energy * energy
        mean_row += probability * row

    assert total == 1
    assert mean_energy == expected_energy
    assert mean_energy2 - mean_energy * mean_energy == expected_variance
    assert mean_row == expected_row


def main() -> None:
    exact_check(5, (Fraction(1, 3), Fraction(-1, 4), Fraction(2, 5), Fraction(1, 7), Fraction(-2, 9)))
    exact_check(6, (Fraction(1, 2), Fraction(1, 5), Fraction(-1, 3), Fraction(2, 7), Fraction(-1, 6), Fraction(3, 8)))
    exact_check(7, (Fraction(-1, 5), Fraction(1, 4), Fraction(2, 9), Fraction(-3, 10), Fraction(1, 7), Fraction(2, 5), Fraction(-1, 8)))
    print("PASS anisotropic_flip_noise_r38_check (exact Fraction identities)")


if __name__ == "__main__":
    main()
