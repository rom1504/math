#!/usr/bin/env python3
"""Verify the exact rational SDP certificate for the Hadamard profile lift."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from phase2b_verify_phi6_collision import matrix_from_code


BASE_CODE = 5850642905
DUAL_NUMERATORS = (445, 490, 661, 668, 436, 645, 405, 427, 485, 513)
DUAL_DENOMINATOR = 200


def bareiss_determinant(matrix: list[list[int]]) -> int:
    """Return an exact determinant using fraction-free elimination."""
    a = [row[:] for row in matrix]
    n = len(a)
    sign = 1
    previous = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            swap = next((i for i in range(k + 1, n) if a[i][k]), None)
            if swap is None:
                return 0
            a[k], a[swap] = a[swap], a[k]
            sign *= -1
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                a[i][j] = (
                    a[i][j] * pivot - a[i][k] * a[k][j]
                ) // previous
        previous = pivot
    return sign * a[-1][-1]


def leading_principal_minors(matrix: np.ndarray) -> list[int]:
    return [
        bareiss_determinant(matrix[:size, :size].tolist())
        for size in range(1, len(matrix) + 1)
    ]


def sylvester_hadamard(order: int) -> np.ndarray:
    matrix = np.ones((1, 1), dtype=np.int64)
    while len(matrix) < order:
        matrix = np.block([[matrix, matrix], [matrix, -matrix]])
    if len(matrix) != order:
        raise ValueError("order must be a power of two")
    return matrix


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    base = matrix_from_code(BASE_CODE)
    q_matrix = base + np.eye(len(base), dtype=np.int64)
    numerators = np.asarray(DUAL_NUMERATORS, dtype=np.int64)

    # Multiplying Diag(y) +/- Q/2 by 200 gives these integer matrices.
    scaled_minus = np.diag(numerators) - (DUAL_DENOMINATOR // 2) * q_matrix
    scaled_plus = np.diag(numerators) + (DUAL_DENOMINATOR // 2) * q_matrix
    minors_minus = leading_principal_minors(scaled_minus)
    minors_plus = leading_principal_minors(scaled_plus)
    if not all(value > 0 for value in minors_minus + minors_plus):
        raise AssertionError("rational SDP dual is not strictly feasible")

    h4 = sylvester_hadamard(4)
    eigenvector = np.asarray((1, 1, 1, -1), dtype=np.int64)
    if not np.array_equal(h4 @ eigenvector, 2 * eigenvector):
        raise AssertionError("H4 Boolean eigenvector certificate failed")
    if int(np.trace(h4)) != 0:
        raise AssertionError("H4 trace certificate failed")

    numerator_sum = int(numerators.sum())
    if numerator_sum >= 26 * DUAL_DENOMINATOR:
        raise AssertionError("dual certificate does not separate 26")

    output = {
        "schema": "quadratic-signing-phase2c-hadamard-sdp-certificate-v1",
        "classification": "exact rational PSD and Boolean-eigenvector certificate",
        "base_code": BASE_CODE,
        "q_definition": "Q=A+I_10",
        "dual_numerators": list(DUAL_NUMERATORS),
        "dual_denominator": DUAL_DENOMINATOR,
        "dual_sum_numerator": numerator_sum,
        "dual_sum": f"{numerator_sum}/{DUAL_DENOMINATOR}",
        "scaled_matrix_definition": (
            "200*(Diag(y) +/- Q/2) = Diag(dual_numerators) +/- 100*Q"
        ),
        "scaled_minus_leading_principal_minors": minors_minus,
        "scaled_plus_leading_principal_minors": minors_plus,
        "sylvester_h4": h4.tolist(),
        "h4_boolean_eigenvector": eigenvector.tolist(),
        "h4_eigenvalue": 2,
        "h4_trace": int(np.trace(h4)),
        "separation": "26-207/8 = 1/8",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps(output, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
