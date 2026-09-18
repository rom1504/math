#!/usr/bin/env python3
"""Exact modular obstruction to inheriting the seed degree-transfer pair.

A full modular column rank certifies full rank over Q(zeta), hence over C.
No floating rank or numerical tolerance enters that implication.
"""
import importlib.util
import json
from pathlib import Path

import numpy as np


def modular_rank(matrix, prime):
    matrix = np.array(matrix, dtype=np.int64, copy=True) % prime
    rank = 0
    for column in range(matrix.shape[1]):
        candidates = np.flatnonzero(matrix[rank:, column])
        if not len(candidates):
            continue
        pivot = rank + int(candidates[0])
        matrix[[rank, pivot]] = matrix[[pivot, rank]]
        matrix[rank] = matrix[rank] * pow(int(matrix[rank, column]), prime - 2, prime) % prime
        factors = matrix[rank + 1:, column].copy()
        matrix[rank + 1:] = (matrix[rank + 1:] - factors[:, None] * matrix[rank][None, :]) % prime
        rank += 1
        if rank == min(matrix.shape):
            break
    return rank


def main():
    source = Path(__file__).with_name("bh_boundedness_counterexamples_2026_09_18_tangent.py")
    spec = importlib.util.spec_from_file_location("tangent", str(source))
    t = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(t)
    roots = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
    masks = [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 16, 17, 18, 20, 24]
    powers = [0, 2, 0, 5, 1, 4, 2, 4, 3, 1, 4, 1, 2, 2, 3, 2]
    f4 = dict(zip(masks, [roots[j] for j in powers]))
    a = {0: (1, -1), 4: (0, 1), 8: (1, 0)}
    fvalues = []
    avalues = []
    for x in range(32):
        fraw = t.evaluate(f4, x)
        fvalues.append((fraw[0] // 4, fraw[1] // 4))
        avalues.append(t.evaluate(a, x))
    # This small prime has zeta=375, a root of z^2-z+1.
    prime = 1009
    assert all(prime % divisor for divisor in range(2, 32))
    zeta = next(z for z in range(prime) if (z*z-z+1) % prime == 0)
    assert zeta == 375 or zeta == 635
    results = []
    for exponent in [10, 104]:
        denominator = 2**exponent
        physical = []
        for y in range(32):
            for x in range(32):
                F = t.mul(fvalues[x], fvalues[y])
                ux = t.mul(avalues[x], t.conj(fvalues[x]))
                uy = t.mul(avalues[y], t.conj(fvalues[y]))
                first = t.mul(ux, t.conj(uy))
                second = t.mul(t.conj(ux), uy)
                quotient = (first[0]-second[0], first[1]-second[1])
                R = 0 if quotient == (0, 0) else 1
                numerator = (denominator*denominator + 48 - 96*R + 2*denominator*quotient[0],
                             2*denominator*quotient[1])
                numerator = t.mul(t.mul(F, F), numerator)
                physical.append(numerator)
        exact_coefficients = list(physical)
        width = 1
        while width < len(exact_coefficients):
            for start in range(0, len(exact_coefficients), 2*width):
                for offset in range(width):
                    i, j = start+offset, start+offset+width
                    left, right = exact_coefficients[i], exact_coefficients[j]
                    exact_coefficients[i] = t.add(left, right)
                    exact_coefficients[j] = (left[0]-right[0], left[1]-right[1])
            width *= 2
        square_degree = max(bin(i).count("1") for i, z in enumerate(exact_coefficients) if z != (0, 0))
        assert square_degree == 8
        coefficients = np.array([(z[0]+zeta*z[1]) % prime for z in exact_coefficients], dtype=np.int64)
        # Common nonzero factors 1024*(denominator^2+48) are irrelevant.
        assert (denominator*denominator + 48) % prime
        weights = [bin(i).count("1") for i in range(1024)]
        for low_degree in range(5):
            high_degree = 8-low_degree
            columns = np.array([i for i in range(1024) if weights[i] <= low_degree])
            rows = np.array([i for i in range(1024) if weights[i] > high_degree])
            matrix = coefficients[rows[:, None] ^ columns[None, :]]
            rank = modular_rank(matrix, prime)
            assert rank == [0, 11, 56, 176, 385][low_degree]
            results.append({"epsilon_exponent": exponent, "low_degree": low_degree,
                            "other_degree": high_degree, "columns": len(columns),
                            "rows": len(rows), "modular_rank": rank,
                            "kernel_dimension_upper_bound": len(columns)-rank})
    print(json.dumps({"status": "PASS", "prime": prime, "zeta": zeta, "square_degree_exact": 8,
                      "reports": results}, indent=2))


if __name__ == "__main__":
    main()
