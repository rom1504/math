#!/usr/bin/env python3
"""Exact finite diagnostics for the Boolean-port Rademacher coreset bound."""

from fractions import Fraction
from itertools import product


def signs(n):
    return list(product((-1, 1), repeat=n))


def projective_signs(n):
    # Fix the first coordinate to +1 as a canonical representative.
    return [(1,) + tail for tail in product((-1, 1), repeat=n - 1)]


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def check_sample(rows, endpoints, sigmas, p):
    absolute_average = Fraction(0)
    linear_average = Fraction(0)
    l1_average = Fraction(0)

    for sigma in sigmas:
        absolute_sup = max(
            abs(sum(sg * abs(dot(row, eps)) for sg, row in zip(sigma, rows)))
            for eps in endpoints
        )
        linear_sup = max(
            abs(sum(sg * dot(row, eps) for sg, row in zip(sigma, rows)))
            for eps in endpoints
        )
        coordinate_l1 = sum(
            abs(sum(sg * row[i] for sg, row in zip(sigma, rows)))
            for i in range(p)
        )
        assert linear_sup == coordinate_l1
        absolute_average += Fraction(absolute_sup, p * len(sigmas))
        linear_average += Fraction(linear_sup, p * len(sigmas))
        l1_average += Fraction(coordinate_l1, p * len(sigmas))

    assert linear_average == l1_average
    assert absolute_average <= 2 * linear_average


def main():
    cases = 0
    for p in range(1, 5):
        rows_space = projective_signs(p)
        endpoints = projective_signs(p)
        for k in range(1, 5):
            sigmas = signs(k)
            for rows in product(rows_space, repeat=k):
                check_sample(rows, endpoints, sigmas, p)
                cases += 1
    print(f"PASS: exact contraction and l1 identities for {cases} sample matrices")


if __name__ == "__main__":
    main()
