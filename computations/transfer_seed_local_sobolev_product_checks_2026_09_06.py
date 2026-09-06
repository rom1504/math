#!/usr/bin/env python3
"""Finite-cube checks of original Sobolev identities and product bounds.

No input files or output artifacts. These checks do not establish an
asymptotic two-root contraction estimate.
"""

from itertools import combinations
import math

import numpy as np


def popcount(value):
    return bin(int(value)).count('1')


def walsh(values):
    result = np.asarray(values, dtype=float).copy()
    size = result.size
    stride = 1
    while stride < size:
        for start in range(0, size, 2 * stride):
            left = result[start:start + stride].copy()
            right = result[start + stride:start + 2 * stride].copy()
            result[start:start + stride] = left + right
            result[start + stride:start + 2 * stride] = left - right
        stride *= 2
    return result


def derivative(values, mask):
    indices = np.arange(values.size)
    signs = np.array([(-1.0) ** popcount(j & mask) for j in indices])
    total = np.zeros_like(values)
    subset = mask
    while True:
        total += (-1) ** popcount(subset) * values[indices ^ subset]
        if subset == 0:
            break
        subset = (subset - 1) & mask
    return signs * total / 2 ** popcount(mask)


def main():
    rng = np.random.default_rng(2026090621)
    sobolev_cases = 0
    product_cases = 0
    coefficient_cases = 0
    for degree in range(1, 13):
        coefficient_constant = math.factorial(degree) * degree ** degree
        for right_degree in range(101):
            for overlap in range(min(degree, right_degree) + 1):
                square = (
                    math.factorial(degree) * math.factorial(right_degree)
                    * math.factorial(degree + right_degree - 2 * overlap)
                    / (math.factorial(overlap) ** 2
                       * math.factorial(degree - overlap) ** 2
                       * math.factorial(right_degree - overlap) ** 2)
                )
                assert square <= coefficient_constant * (right_degree + 1) ** degree * (1 + 1e-12)
                coefficient_cases += 1
    for n in range(3, 9):
        size = 1 << n
        degrees = np.array([popcount(mask) for mask in range(size)])
        for _ in range(5):
            source = rng.normal(size=size) * (degrees <= 3)
            source /= np.linalg.norm(source)
            w = walsh(source)
            a = np.cos(w) + 0.3 * np.sin(2 * w)
            ahat = walsh(a) / size
            for order in range(1, 4):
                spectral = sum(math.comb(int(q), order) * c * c
                               for q, c in zip(degrees, ahat) if q >= order)
                literal = 0.0
                for slots in combinations(range(n), order):
                    mask = sum(1 << j for j in slots)
                    literal += float(np.mean(derivative(a, mask) ** 2))
                assert abs(spectral - literal) < 1e-11
                sobolev_cases += 1
            for degree in range(1, 4):
                xhat = rng.normal(size=size) * (degrees == degree)
                xhat /= np.linalg.norm(xhat)
                x = walsh(xhat)
                yhat = rng.normal(size=size)
                y = walsh(yhat)
                constant = (degree + 1) ** 2 * math.factorial(degree) * degree ** degree
                lhs = float(np.mean((x * y) ** 2))
                rhs = constant * float(np.sum((degrees + 1) ** degree * yhat ** 2))
                assert lhs <= rhs * (1 + 1e-12)
                product_cases += 1
    print('exact_sobolev_identity_cases', sobolev_cases)
    print('finite_cube_weighted_product_cases', product_cases)
    print('factorial_coefficient_cases', coefficient_cases)
    print('PASS: finite local algebra; no open two-root conclusion')


if __name__ == '__main__':
    main()
