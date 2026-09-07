#!/usr/bin/env python3
"""Independent exact audit of the balanced-partition second moment."""

from fractions import Fraction
from itertools import combinations, permutations, product


def energy_parts(A, z, parts):
    return sum(
        2 * A[i][j] * z[i] * z[j]
        for block in parts
        for i, j in combinations(block, 2)
    )


def equipartitions(vertices, s):
    vertices = tuple(vertices)
    if not vertices:
        yield ()
        return
    first = vertices[0]
    for tail in combinations(vertices[1:], s - 1):
        block = (first,) + tail
        rest = tuple(v for v in vertices if v not in block)
        for suffix in equipartitions(rest, s):
            yield (block,) + suffix


def formula(A, z, s):
    n = len(A)
    p = Fraction(s - 1, n - 1)
    ca = -p * Fraction(n - s, (n - 1) * (n - 2))
    cd = 2 * p * Fraction(n - s, (n - 1) * (n - 2) * (n - 3))
    h = sum(A[i][j] * z[i] * z[j] for i in range(n) for j in range(n))
    a2 = sum(
        sum(A[i][j] * z[j] for j in range(n)) ** 2
        for i in range(n)
    )
    # Unordered-edge weights are w_ij=2 A_ij z_i z_j.
    first = 2 * n * (n - 1) * p * (1 - p)
    adjacent = 2 * a2 - 2 * n * (n - 1)
    disjoint = Fraction(h * h, 2) - 2 * a2 + n * (n - 1)
    return p * h, first + 2 * ca * adjacent + 2 * cd * disjoint


def audit(n, s):
    parts = tuple(equipartitions(range(n), s))
    edges = tuple(combinations(range(n), 2))
    # A few deterministic signing patterns, including all-positive.
    patterns = [
        tuple(1 for _ in edges),
        tuple(1 if k % 2 else -1 for k, _ in enumerate(edges)),
        tuple(1 if (i + j) % 3 else -1 for i, j in edges),
    ]
    for signs in patterns:
        A = [[0] * n for _ in range(n)]
        for (i, j), a in zip(edges, signs):
            A[i][j] = A[j][i] = a
        for z in product((-1, 1), repeat=n):
            vals = tuple(energy_parts(A, z, P) for P in parts)
            mean = Fraction(sum(vals), len(vals))
            var = sum((Fraction(v) - mean) ** 2 for v in vals) / len(vals)
            expected_mean, expected_var = formula(A, z, s)
            assert mean == expected_mean
            assert var == expected_var
            assert var <= 4 * n * s
    return len(parts)


if __name__ == "__main__":
    assert audit(6, 2) == 15
    assert audit(6, 3) == 10
    assert audit(8, 2) == 105
    print("PASS: exact equipartition mean/variance and sampled Var <= 4ns")
