#!/usr/bin/env python3
"""Exact rational checks for the Wave 38 complement-flip noise identities."""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product


def matvec(a: list[list[int]], x: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sum(a[i][j] * x[j] for j in range(len(x))) for i in range(len(x)))


def quad(a: list[list[int]], x: tuple[int, ...]) -> int:
    ax = matvec(a, x)
    return sum(x[i] * ax[i] for i in range(len(x)))


def row_square(a: list[list[int]], x: tuple[int, ...]) -> int:
    return sum(v * v for v in matvec(a, x))


def complement_matrix(a: list[list[int]], selector: frozenset[int]) -> list[list[int]]:
    n = len(a)
    return [
        [
            0 if i == j else a[i][j] * (1 if i in selector and j in selector else -1)
            for j in range(n)
        ]
        for i in range(n)
    ]


def deterministic_signing(n: int) -> list[list[int]]:
    # Symmetric, zero diagonal, and deliberately nonsymmetric in its row sums.
    a = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            value = 1 if ((i * j + i + 2 * j + j * j) % 5) < 3 else -1
            a[i][j] = a[j][i] = value
    return a


def exact_product_check(n: int, flip_probability: Fraction) -> None:
    a_mat = deterministic_signing(n)
    selector = frozenset(range(1, n))
    b_mat = complement_matrix(a_mat, selector)

    spins = list(product((-1, 1), repeat=n))
    ground = max(spins, key=lambda x: abs(quad(b_mat, x)))
    raw = quad(b_mat, ground)
    sigma = 1 if raw >= 0 else -1
    m_value = sigma * raw

    # C = sigma D_x B D_x has all-ones ground energy M.
    c_mat = [
        [sigma * ground[i] * b_mat[i][j] * ground[j] for j in range(n)]
        for i in range(n)
    ]
    r_b = sum(sum(row) ** 2 for row in c_mat)
    r_a = row_square(a_mat, ground)

    eps = flip_probability
    rho = 1 - 2 * eps
    alpha = rho * rho

    mean_e = Fraction(0)
    mean_e2 = Fraction(0)
    mean_r = Fraction(0)
    total = Fraction(0)
    for bits in product((0, 1), repeat=n):
        k = sum(bits)
        probability = eps**k * (1 - eps) ** (n - k)
        z = tuple(-1 if bit else 1 for bit in bits)
        noisy = tuple(ground[i] * z[i] for i in range(n))
        energy = sigma * quad(b_mat, noisy)
        row = row_square(a_mat, noisy)
        total += probability
        mean_e += probability * energy
        mean_e2 += probability * energy * energy
        mean_r += probability * row

    assert total == 1
    expected_mean = alpha * m_value
    expected_variance = (
        4 * alpha * (1 - alpha) * r_b
        + 2 * (1 - alpha) ** 2 * n * (n - 1)
    )
    expected_row = alpha * r_a + (1 - alpha) * n * (n - 1)
    assert mean_e == expected_mean
    assert mean_e2 - mean_e * mean_e == expected_variance
    assert mean_r == expected_row

    # Conditional on an exact Hamming radius, both quadratic expectations
    # have the same pair-correlation coefficient alpha_k.
    for k in range(n + 1):
        energies: list[int] = []
        rows: list[int] = []
        for flipped in combinations(range(n), k):
            flipped_set = set(flipped)
            z = tuple(-1 if i in flipped_set else 1 for i in range(n))
            noisy = tuple(ground[i] * z[i] for i in range(n))
            energies.append(sigma * quad(b_mat, noisy))
            rows.append(row_square(a_mat, noisy))
        alpha_k = Fraction((n - 2 * k) ** 2 - n, n * (n - 1))
        assert Fraction(sum(energies), len(energies)) == alpha_k * m_value
        assert Fraction(sum(rows), len(rows)) == (
            alpha_k * r_a + (1 - alpha_k) * n * (n - 1)
        )


def main() -> None:
    for n, eps in ((5, Fraction(1, 3)), (6, Fraction(2, 5)), (7, Fraction(3, 7))):
        exact_product_check(n, eps)
    print("PASS dual_flip_compression_r38_check (all identities exact over Fraction)")


if __name__ == "__main__":
    main()
