#!/usr/bin/env python3
"""Independent numerical smoke test of the 2026-09-05 polar-Gram bounds.

The proofs are in artifacts/fresh_limit_variational_2026_09_05.md.
SDP values here are numerical diagnostics, never exact certificates.
"""

from __future__ import annotations

import json
import math

import cvxpy as cp
import numpy as np

from conference_prime_square import PrimeSquare


def exact_cap(matrix: np.ndarray) -> int:
    n = len(matrix)
    codes = np.arange(1 << (n - 1), dtype=np.uint64)
    spins = np.ones((len(codes), n), dtype=np.int64)
    shifts = np.arange(n - 1, dtype=np.uint64)
    spins[:, 1:] -= 2 * ((codes[:, None] >> shifts) & 1).astype(np.int64)
    energy = np.einsum("bi,bi->b", spins @ matrix, spins) // 2
    return int(np.abs(energy).max())


def paley_prime(prime: int) -> np.ndarray:
    assert prime % 4 == 1
    squares = {i * i % prime for i in range(1, prime)}
    matrix = np.zeros((prime + 1, prime + 1), dtype=np.int64)
    matrix[0, 1:] = matrix[1:, 0] = 1
    for i in range(prime):
        for j in range(prime):
            if i != j:
                matrix[i + 1, j + 1] = 1 if (i - j) % prime in squares else -1
    assert np.array_equal(matrix @ matrix, prime * np.eye(prime + 1))
    return matrix


def diagnostic(matrix: np.ndarray, label: str) -> dict[str, object]:
    n = len(matrix)
    values, vectors = np.linalg.eigh(matrix)
    absolute = (vectors * np.abs(values)) @ vectors.T
    h = np.diag(absolute)
    nuclear = float(np.sum(np.abs(values)))
    denominator = np.sqrt(h[:, None] * h[None, :])
    x_plus = (absolute + matrix) / denominator
    x_minus = (absolute - matrix) / denominator
    polar_vector_bound = float(np.sum(matrix * matrix / denominator) / 2)
    rounded_plus = float(np.sum(matrix * np.arcsin(np.clip(x_plus, -1, 1))) / math.pi)
    rounded_minus = float(np.sum(matrix * np.arcsin(np.clip(x_minus, -1, 1))) / math.pi)
    paired_rounding = (rounded_plus - rounded_minus) / 2
    nuclear_rounding = n * (n - 1) / math.pi * math.asin(min(1.0, n / nuclear))
    cap = exact_cap(matrix)
    variable = cp.Variable((n, n), symmetric=True)
    constraints = [variable >> 0, cp.diag(variable) == 1]
    sdp_values = []
    for sign in (1, -1):
        problem = cp.Problem(cp.Maximize(sign * cp.sum(cp.multiply(matrix, variable)) / 2), constraints)
        sdp_values.append(float(problem.solve(solver="CLARABEL")))
    vector_value = max(sdp_values)
    base = n * math.sqrt(n - 1) / 2
    epsilon = max(0.0, vector_value / base - 1)
    spectral_defect = float(np.mean((np.abs(values) / math.sqrt(n - 1) - 1) ** 2))
    tolerance = 1e-5 * n * n
    assert np.min(np.linalg.eigvalsh(x_plus)) >= -1e-9
    assert np.min(np.linalg.eigvalsh(x_minus)) >= -1e-9
    assert cap + tolerance >= paired_rounding >= nuclear_rounding - tolerance
    assert vector_value + tolerance >= polar_vector_bound
    assert vector_value * nuclear + tolerance >= n * n * (n - 1) / 2
    assert spectral_defect <= 2 * epsilon / (1 + epsilon) + tolerance
    return {
        "label": label,
        "n": n,
        "exact_boolean_cap": cap,
        "numerical_vector_value": vector_value,
        "normalized_nuclear_mass": nuclear / (n * math.sqrt(n - 1)),
        "spectral_defect": spectral_defect,
        "paired_gaussian_bound": paired_rounding,
        "nuclear_arcsine_bound": nuclear_rounding,
        "boolean_vector_ratio": cap / vector_value,
    }


def main() -> None:
    generator = np.random.default_rng(20260905)
    records = []
    for n in range(2, 13):
        for sample in range(3):
            upper = np.triu(generator.choice((-1, 1), (n, n)), 1)
            records.append(diagnostic(upper + upper.T, f"random_{n}_{sample}"))
    for prime in (5, 13, 17):
        records.append(diagnostic(paley_prime(prime), f"conference_prime_{prime}"))
    records.append(diagnostic(PrimeSquare(3).conference().astype(np.int64), "conference_q9"))
    print(json.dumps({"seed": 20260905, "records": records}, indent=2))


if __name__ == "__main__":
    main()
