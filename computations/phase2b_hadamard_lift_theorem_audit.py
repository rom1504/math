#!/usr/bin/env python3
"""Independent exact audit of the equal-profile Hadamard-lift theorem."""

from __future__ import annotations

import argparse
import itertools
import json
from fractions import Fraction
from pathlib import Path

import numpy as np


CODES = (5850642905, 28771662001)
DUAL_NUMERATORS = (445, 490, 661, 668, 436, 645, 405, 427, 485, 513)
REFINED_DUAL_DENOMINATOR = 100_000
REFINED_DUAL_NUMERATORS = (
    222237, 244660, 330058, 333605, 217787,
    322119, 202501, 213301, 242265, 256161,
)


def matrix_from_code(code: int, order: int = 10) -> np.ndarray:
    matrix = np.ones((order, order), dtype=np.int64)
    np.fill_diagonal(matrix, 0)
    bit = 0
    for i in range(1, order):
        for j in range(i + 1, order):
            if code & (1 << bit):
                matrix[i, j] = matrix[j, i] = -1
            bit += 1
    return matrix


def exact_ldl_pivots(matrix: list[list[int]]) -> list[Fraction]:
    """Return exact diagonal pivots in an unpivoted LDL^T factorization."""
    order = len(matrix)
    lower = [[Fraction(0) for _ in range(order)] for _ in range(order)]
    diagonal: list[Fraction] = []
    for j in range(order):
        pivot = Fraction(matrix[j][j]) - sum(
            lower[j][k] * lower[j][k] * diagonal[k] for k in range(j)
        )
        if not pivot:
            raise AssertionError("zero LDL pivot")
        diagonal.append(pivot)
        lower[j][j] = Fraction(1)
        for i in range(j + 1, order):
            lower[i][j] = (
                Fraction(matrix[i][j])
                - sum(
                    lower[i][k] * lower[j][k] * diagonal[k]
                    for k in range(j)
                )
            ) / pivot
    return diagonal


def sylvester(order: int) -> np.ndarray:
    matrix = np.ones((1, 1), dtype=np.int64)
    while len(matrix) < order:
        matrix = np.block([[matrix, matrix], [matrix, -matrix]])
    return matrix


def positive_base_maximum(matrix: np.ndarray) -> tuple[int, list[int]]:
    best = -10**9
    witness: list[int] = []
    for tail in itertools.product((-1, 1), repeat=len(matrix) - 1):
        spin = np.asarray((1, *tail), dtype=np.int64)
        energy = int(spin @ matrix @ spin // 2)
        if energy > best:
            best = energy
            witness = spin.tolist()
    return best, witness


def lift(base: np.ndarray, hadamard: np.ndarray) -> np.ndarray:
    diagonal = np.diag(np.diag(hadamard))
    return np.kron(base, hadamard) + np.kron(
        np.eye(len(base), dtype=np.int64), hadamard - diagonal
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    bases = [matrix_from_code(code) for code in CODES]
    first_q = bases[0] + np.eye(10, dtype=np.int64)
    certificates = {}
    for sign in (-1, 1):
        integer_matrix = (
            np.diag(np.asarray(DUAL_NUMERATORS, dtype=np.int64))
            + sign * 100 * first_q
        )
        pivots = exact_ldl_pivots(integer_matrix.tolist())
        if not all(pivot > 0 for pivot in pivots):
            raise AssertionError("dual matrix is not positive definite")
        certificates[str(sign)] = [str(pivot) for pivot in pivots]

    refined_certificates = {}
    for sign in (-1, 1):
        # 2D*(Diag(y) + sign*Q/2), for one common refined y.
        integer_matrix = (
            2 * np.diag(np.asarray(REFINED_DUAL_NUMERATORS, dtype=np.int64))
            + sign * REFINED_DUAL_DENOMINATOR * first_q
        )
        pivots = exact_ldl_pivots(integer_matrix.tolist())
        if not all(pivot > 0 for pivot in pivots):
            raise AssertionError("refined common dual matrix is not positive definite")
        refined_certificates[str(sign)] = [str(pivot) for pivot in pivots]

    h4 = sylvester(4)
    seed = np.asarray((-1, -1, -1, 1), dtype=np.int64)
    if not np.array_equal(h4 @ seed, 2 * seed):
        raise AssertionError("Boolean positive eigenvector failed")
    if int(np.trace(h4)) != 0 or not np.array_equal(h4 @ h4, 4 * np.eye(4)):
        raise AssertionError("Hadamard identities failed")

    positive_maxima = []
    lift_witnesses = []
    for base in bases:
        maximum, sigma = positive_base_maximum(base)
        positive_maxima.append(maximum)
        signing = lift(base, h4)
        off_diagonal = signing[~np.eye(len(signing), dtype=bool)]
        if np.any(np.diag(signing)) or not np.all(np.abs(off_diagonal) == 1):
            raise AssertionError("invalid lifted signing")
        spin = np.kron(np.asarray(sigma, dtype=np.int64), seed)
        energy = int(spin @ signing @ spin // 2)
        predicted = (maximum + 5) * 8
        if energy != predicted:
            raise AssertionError((energy, predicted))
        lift_witnesses.append({
            "base_positive_maximum": maximum,
            "base_spin": sigma,
            "lift_energy_at_k4": energy,
        })

    coefficient = Fraction(sum(DUAL_NUMERATORS), 200)
    if coefficient != Fraction(207, 8):
        raise AssertionError(coefficient)
    gap = Fraction(26) - coefficient
    refined_coefficient = Fraction(
        sum(REFINED_DUAL_NUMERATORS), REFINED_DUAL_DENOMINATOR
    )
    refined_gap = Fraction(26) - refined_coefficient
    output = {
        "schema": "quadratic-signing-phase2b-hadamard-theorem-audit-v1",
        "classification": "independent exact rational and tensor-factor audit",
        "dual_interpretation": (
            "y_i=dual_numerator_i/200; exact LDL proves "
            "Diag(y) plus/minus (A+I)/2 positive definite"
        ),
        "exact_ldl_pivots_after_scaling_by_200": certificates,
        "dual_objective": str(coefficient),
        "refined_common_dual": {
            "denominator": REFINED_DUAL_DENOMINATOR,
            "numerators": list(REFINED_DUAL_NUMERATORS),
            "objective": str(refined_coefficient),
            "gap_coefficient": str(refined_gap),
            "exact_ldl_pivots_for_both_signs": refined_certificates,
        },
        "base_positive_maxima": positive_maxima,
        "h4_trace": int(np.trace(h4)),
        "h4_positive_boolean_eigenvector": seed.tolist(),
        "k4_lift_witnesses": lift_witnesses,
        "first_lift_cap_upper_coefficient": str(coefficient),
        "second_lift_cap_lower_coefficient": "26",
        "gap_coefficient": str(gap),
        "normalized_order_N_gap": "1/(8*10^(3/2))",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps(output, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
