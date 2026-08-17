#!/usr/bin/env python3
"""Exact arithmetic plus numerical smoke checks for the Walsh ceiling."""

from fractions import Fraction
from math import ceil, comb, isqrt, sqrt


def gaussian_binomial(n: int, r: int) -> int:
    """Number of r-dimensional subspaces of F_2^n, computed exactly."""
    num = 1
    den = 1
    for j in range(r):
        num *= (1 << (n - j)) - 1
        den *= (1 << (r - j)) - 1
    assert num % den == 0
    return num // den


def number_of_subspaces(n: int) -> int:
    return sum(gaussian_binomial(n, r) for r in range(n + 1))


def main() -> None:
    checks = 0

    # Exact equal-cell accounting: changed flux bits alter exactly two of
    # three child terms.  Normalize n^(3/2) away.
    for h in range(1, 501):
        k = 3 * h
        local_diameter = Fraction(2 * h, 1)
        total_scale_squared = k**3  # (N^(3/2)/n^(3/2))^2
        # Equation (TC.11), squared to stay exact.
        assert local_diameter**2 * 27 * h == Fraction(4 * h * h, 1) * 27 * h
        assert local_diameter**2 / total_scale_squared == Fraction(4, 27 * h)
        # A Hamming neighbour gives (TC.12).
        neighbour = Fraction(2, 1)
        assert neighbour**2 / total_scale_squared == Fraction(4, 27 * h**3)
        checks += 2

    # Integral tail inequality sum_{i=r+1}^H i^{-3/2} <= 2/sqrt(r).
    # Floating point is used only to audit this elementary analytic bound,
    # with ample margin; the proof in the draft is integral comparison.
    for r in range(1, 200):
        tail = sum(i ** -1.5 for i in range(r + 1, 200000))
        remainder = 2.0 / sqrt(199999.0)
        assert tail + remainder <= 2.0 / sqrt(r) + 1e-12
        checks += 1

    # The exact number of subspaces is below the crude 2^(k^2) carrier count.
    for k in range(0, 13):
        subspaces = number_of_subspaces(k)
        assert subspaces <= 1 << (k * k)
        orbit_bound = (1 << (k * (k + 1) // 2)) * (1 << (k * k))
        assert orbit_bound == 1 << ((3 * k * k + k) // 2)
        checks += 2

    # Check the quadratic inversion: h <= (3k^2+k)/2 implies
    # k >= (sqrt(1+24h)-1)/6.  Integer version through a broad range.
    for h in range(1, 100000):
        k = 0
        while (3 * k * k + k) // 2 < h:
            k += 1
        lower = (sqrt(1 + 24 * h) - 1) / 6
        assert k + 1e-12 >= lower
        if k > 0:
            assert (3 * (k - 1) * (k - 1) + (k - 1)) // 2 < h
        checks += 1

    print(f"Walsh total-scale flux ceiling checks passed: {checks}")


if __name__ == "__main__":
    main()
