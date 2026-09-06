#!/usr/bin/env python3
"""Exact finite-difference cover and numerical matrix-bound checks.

These test finite algebra, not any asymptotic local feedback theorem.
No input files or generated outputs; run with .venv/bin/python.
"""

from itertools import combinations, product
import math

import numpy as np


def op(a):
    return float(np.linalg.norm(a, 2))


def subsets(items):
    return [tuple(items[j] for j in range(len(items)) if bits >> j & 1)
            for bits in range(1 << len(items))]


def matrix_checks(rng):
    count = 0
    for rows in range(1, 9):
        for cols in range(1, 9):
            for scale in (0.1, 1, 10):
                j = scale * rng.normal(size=(rows, cols))
                t = rng.uniform(-3, 3)
                lhs = op(np.exp(1j*t*j)-1)
                rhs = abs(t)*op(j)+(t*t/2)*op(j)**2
                assert lhs <= rhs + 1e-10*max(1, rhs)
                count += 1
    return count


def cover_checks(rng):
    n = 5
    rows = 4
    all_sets = subsets(tuple(range(n)))
    coeff = {u: rng.normal(size=rows)/(1+len(u)) for u in all_sets if len(u) <= 3}

    def value(s):
        return sum(c * np.prod(s[list(u)]) for u, c in coeff.items())

    def derivative(s, u):
        uset = set(u)
        return sum((c * np.prod(s[list(set(v)-uset)]) for v, c in coeff.items()
                    if uset <= set(v)), start=np.zeros(rows))

    t = 0.73
    count = 0
    operator_count = 0
    for signs in product((-1.0, 1.0), repeat=n):
        s = np.array(signs)
        w = value(s)
        derivatives = {}
        a_matrices = {}
        for degree in (1, 2, 3):
            columns = list(combinations(range(n), degree))
            j = np.stack([derivative(s, u) for u in columns], axis=1)
            derivatives[degree] = j
            chars = np.array([np.prod(s[list(u)]) for u in columns])
            a_matrices[degree] = np.exp(1j*t*(-2)**degree*j*chars)-1
            expected_bound = 2**degree*abs(t)*op(j)+2**(2*degree-1)*t*t*op(j)**2
            assert op(a_matrices[degree]) <= expected_bound+1e-9*max(1, expected_bound)
            operator_count += 1

        for degree in (1, 2, 3):
            for u in combinations(range(n), degree):
                parts = subsets(u)
                direct = np.zeros(rows, dtype=complex)
                for v in parts:
                    flipped = s.copy()
                    flipped[list(v)] *= -1
                    direct += (-1)**len(v)*np.exp(1j*t*value(flipped))
                direct *= np.prod(s[list(u)])/2**degree

                nonempty = parts[1:]
                a_values = [np.exp(1j*t*(-2)**len(v)*np.prod(s[list(v)])*derivative(s, v))-1
                            for v in nonempty]
                cover = np.zeros(rows, dtype=complex)
                for bits in range(1, 1 << len(nonempty)):
                    family = [j for j in range(len(nonempty)) if bits >> j & 1]
                    union = set().union(*(set(nonempty[j]) for j in family))
                    if union == set(u):
                        cover += np.prod([a_values[j] for j in family], axis=0)
                cover *= (-1)**degree*np.prod(s[list(u)])*np.exp(1j*t*w)/2**degree
                assert np.max(np.abs(direct-cover)) < 2e-12
                count += rows
    return count, operator_count


def main():
    rng = np.random.default_rng(2026090607)
    print('entrywise_exponential_matrix_cases', matrix_checks(rng))
    cases, operators = cover_checks(rng)
    print('exact_boolean_cover_scalar_cases', cases)
    print('modified_derivative_matrix_cases', operators)
    print('PASS: finite algebra only; no small proper BC cut asserted')


if __name__ == '__main__':
    main()
