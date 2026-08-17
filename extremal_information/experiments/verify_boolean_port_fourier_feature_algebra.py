#!/usr/bin/env python3
"""Exact checks for the Boolean-port Fourier feature algebra."""

from __future__ import annotations

from fractions import Fraction
from itertools import product
from math import comb, prod

import numpy as np


def double_factorial(k: int) -> int:
    if k <= 0:
        return 1
    return prod(range(k, 0, -2))


def predicted_multiplier(p: int, degree: int) -> Fraction:
    P = p if p % 2 == 0 else p + 1
    if degree == 0:
        return Fraction(double_factorial(P - 1), double_factorial(P - 2))
    k = degree // 2
    return Fraction(
        (-1) ** (k - 1)
        * double_factorial(2 * k - 3)
        * double_factorial(P - 2 * k - 1),
        double_factorial(P - 2),
    )


def exact_multiplier(p: int, degree: int) -> Fraction:
    total = 0
    for z in product((-1, 1), repeat=p):
        total += abs(sum(z)) * prod(z[:degree])
    return Fraction(total, 2**p)


def projective_rows(p: int) -> np.ndarray:
    # Gauge the first coordinate to +1.
    return np.asarray([(1, *tail) for tail in product((-1, 1), repeat=p - 1)], dtype=int)


def response_matrix(p: int) -> np.ndarray:
    rows = projective_rows(p)
    return np.abs(rows @ rows.T)


def projective_character_matrix(p: int) -> tuple[np.ndarray, list[int]]:
    """Exact character table, with each character's even full-cube degree."""
    rows = projective_rows(p)
    subsets = list(product((0, 1), repeat=p - 1))
    table = []
    degrees = []
    for subset in subsets:
        tail_degree = sum(subset)
        # If the tail monomial has odd degree, multiply by coordinate 0 to
        # obtain its unique even representative on the full cube.  The
        # gauge z_0=+1 makes the two characters identical on G_p.
        degrees.append(tail_degree + (tail_degree & 1))
        row = []
        for z in rows:
            value = 1
            for i, used in enumerate(subset, start=1):
                if used:
                    value *= int(z[i])
            row.append(value)
        table.append(row)
    return np.asarray(table, dtype=np.int64), degrees


def weak_compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for tail in weak_compositions(total - first, parts - 1):
            yield (first, *tail)


def main() -> None:
    checks = 0
    for p in range(1, 13):
        for degree in range(0, p + 1, 2):
            got = exact_multiplier(p, degree)
            want = predicted_multiplier(p, degree)
            assert got == want and got != 0, (p, degree, got, want)
            checks += 1

    # The convolution operator is invertible on the projective cube.  Check
    # this by exact character diagonalization, not a floating determinant.
    for p in range(1, 8):
        R = response_matrix(p)
        C, degrees = projective_character_matrix(p)
        group_size = 2 ** (p - 1)
        assert np.array_equal(C @ C.T, group_size * np.eye(group_size, dtype=np.int64))
        diagonalized = C @ R @ C.T
        assert np.array_equal(diagonalized, np.diag(np.diag(diagonalized)))
        for entry, degree in zip(np.diag(diagonalized), degrees):
            assert Fraction(int(entry), group_size**2) == predicted_multiplier(p, degree)
            assert int(entry) != 0
        checks += group_size + 2

        # PF.2 point-mass response distance is exactly twice projective
        # Hamming distance, after normalization by p.
        rows = projective_rows(p)
        for i in range(group_size):
            for j in range(i + 1, group_size):
                raw_hamming = int(np.count_nonzero(rows[i] != rows[j]))
                projective_hamming = min(raw_hamming, p - raw_hamming)
                response_distance = int(np.max(np.abs(R[i] - R[j])))
                assert Fraction(response_distance, p) == Fraction(
                    2 * projective_hamming, p
                )
                checks += 1

    # At small sizes, every histogram has a distinct labelled response.
    for p, n in ((2, 4), (3, 3), (4, 2)):
        R = response_matrix(p)
        seen: set[tuple[int, ...]] = set()
        count = 0
        for mu in weak_compositions(n, len(R)):
            response = tuple((np.asarray(mu, dtype=int) @ R).tolist())
            assert response not in seen
            seen.add(response)
            count += 1
        assert count == comb(n + 2 ** (p - 1) - 1, 2 ** (p - 1) - 1)
        checks += count + 1

    # Uniform versus even-parity four-port systems.
    rows = np.asarray(list(product((-1, 1), repeat=4)), dtype=int)
    even = rows[np.prod(rows, axis=1) == 1]
    U = rows
    V = np.repeat(even, 2, axis=0)
    eps = rows
    assert np.array_equal(U.T @ U, V.T @ V)
    assert int(np.max(np.sum(np.abs(eps @ U.T), axis=1))) == 24
    assert int(np.max(np.sum(np.abs(eps @ V.T), axis=1))) == 32
    checks += 3

    print(f"Boolean-port Fourier feature checks passed: {checks}")


if __name__ == "__main__":
    main()
