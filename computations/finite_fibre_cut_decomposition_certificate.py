#!/usr/bin/env python3
"""Audit a recurrence-first finite-fibre cut-decomposition certificate.

For a symmetric full-sign microkernel R of order k, index the antipodal
Boolean microtypes u by fixing u_0=1 and set K[u,v]=u^T R v.  The certificate
constant is

    L(R) = min sum_h |c_h|  subject to
           K = sum_h c_h h h^T,

where h ranges over {-1,0,1}^{2^(k-1)} modulo h -> -h.  A certificate gives
cap(A tensor R) <= L(R) cap(A) for every zero-diagonal signing A.

The script exhausts every R for k=2,3,4, solves the finite LP, rationalizes
and independently checks the best primal/dual certificates, and verifies the
closed-form range obstruction used in the accompanying report.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


def antipodal_vectors(k: int) -> np.ndarray:
    return np.asarray(
        [(1,) + tail for tail in itertools.product((-1, 1), repeat=k - 1)],
        dtype=np.int64,
    )


def ternary_atoms(q: int) -> np.ndarray:
    rows = []
    for row in itertools.product((-1, 0, 1), repeat=q):
        if not any(row):
            continue
        if next(value for value in row if value) > 0:
            rows.append(row)
    return np.asarray(rows, dtype=np.int8)


def symmetric_pairs(q: int) -> list[tuple[int, int]]:
    return [(i, j) for i in range(q) for j in range(i, q)]


def atom_matrix(atoms: np.ndarray, pairs: list[tuple[int, int]]) -> np.ndarray:
    return np.asarray(
        [[int(h[i]) * int(h[j]) for h in atoms] for i, j in pairs],
        dtype=np.float64,
    )


def full_sign_matrices(k: int):
    count = k * (k + 1) // 2
    for bits in itertools.product((-1, 1), repeat=count):
        matrix = np.zeros((k, k), dtype=np.int64)
        cursor = 0
        for i in range(k):
            for j in range(i, k):
                matrix[i, j] = matrix[j, i] = bits[cursor]
                cursor += 1
        yield matrix


def upper_triangle(matrix: np.ndarray) -> list[int]:
    return [
        int(matrix[i, j])
        for i in range(len(matrix))
        for j in range(i, len(matrix))
    ]


def solve_lp(
    kernel: np.ndarray,
    atom_columns: np.ndarray,
    pairs: list[tuple[int, int]],
) -> tuple[float, np.ndarray, np.ndarray]:
    target = np.asarray([kernel[i, j] for i, j in pairs], dtype=np.float64)
    equality = np.concatenate((atom_columns, -atom_columns), axis=1)
    result = linprog(
        np.ones(equality.shape[1]),
        A_eq=equality,
        b_eq=target,
        bounds=(0, None),
        method="highs",
    )
    if not result.success:
        raise AssertionError(result.message)
    atom_count = atom_columns.shape[1]
    signed = result.x[:atom_count] - result.x[atom_count:]
    return float(result.fun), signed, result.eqlin.marginals


def rationalize(values: np.ndarray, denominator: int = 10000) -> list[Fraction]:
    return [Fraction(float(value)).limit_denominator(denominator) for value in values]


def verify_certificate(
    kernel: np.ndarray,
    atoms: np.ndarray,
    pairs: list[tuple[int, int]],
    primal_float: np.ndarray,
    dual_float: np.ndarray,
) -> dict[str, object]:
    primal = rationalize(primal_float)
    dual = rationalize(dual_float)
    atom_columns = [
        [int(h[i]) * int(h[j]) for h in atoms] for i, j in pairs
    ]
    target = [int(kernel[i, j]) for i, j in pairs]

    for edge, value in enumerate(target):
        reconstructed = sum(
            primal[column] * atom_columns[edge][column]
            for column in range(len(atoms))
        )
        if reconstructed != value:
            raise AssertionError((edge, reconstructed, value))

    dual_atom_values = []
    for column in range(len(atoms)):
        pairing = sum(
            dual[edge] * atom_columns[edge][column]
            for edge in range(len(pairs))
        )
        if abs(pairing) > 1:
            raise AssertionError((column, pairing))
        dual_atom_values.append(pairing)

    primal_objective = sum(abs(value) for value in primal)
    dual_objective = sum(value * weight for value, weight in zip(target, dual))
    if primal_objective != dual_objective:
        raise AssertionError((primal_objective, dual_objective))

    return {
        "objective": str(primal_objective),
        "primal_support": [
            {
                "atom_index": index,
                "coefficient": str(value),
                "atom": [int(entry) for entry in atoms[index]],
            }
            for index, value in enumerate(primal)
            if value
        ],
        "dual_support": [
            {
                "pair": list(pairs[index]),
                "coefficient": str(value),
            }
            for index, value in enumerate(dual)
            if value
        ],
        "maximum_absolute_dual_atom_pairing": str(
            max(abs(value) for value in dual_atom_values)
        ),
    }


def exact_rademacher_absolute_mean(s: int) -> Fraction:
    return sum(
        Fraction(math.comb(s, j) * abs(s - 2 * j), 2**s)
        for j in range(s + 1)
    )


def range_lower_bound(k: int) -> dict[str, object]:
    # Choosing the larger shore r optimizes the elementary random-y bound.
    candidates = []
    for r in range(1, k):
        s = k - r
        mean = exact_rademacher_absolute_mean(s)
        candidates.append((4 * r * mean, r, s, mean))
    bound, r, s, mean = max(candidates)
    return {
        "split": [r, s],
        "exact_rademacher_absolute_mean": str(mean),
        "exact_energy_range_lower_bound": str(bound),
        "strictly_exceeds_k_to_three_halves": bool(float(bound) > k**1.5),
        "ratio_to_k_to_three_halves": float(bound) / k**1.5,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    orders = []
    expected_minima = {2: 4, 3: 9, 4: 16}
    for k in (2, 3, 4):
        types = antipodal_vectors(k)
        atoms = ternary_atoms(len(types))
        pairs = symmetric_pairs(len(types))
        columns = atom_matrix(atoms, pairs)
        histogram: dict[str, int] = {}
        best_value = float("inf")
        best_data = None
        matrix_count = 0
        for matrix in full_sign_matrices(k):
            kernel = types @ matrix @ types.T
            value, primal, dual = solve_lp(kernel, columns, pairs)
            rounded = round(value, 8)
            histogram[str(rounded)] = histogram.get(str(rounded), 0) + 1
            matrix_count += 1
            if value < best_value - 1e-7:
                best_value = value
                best_data = (matrix.copy(), kernel.copy(), primal, dual)

        if round(best_value) != expected_minima[k] or best_data is None:
            raise AssertionError((k, best_value))
        matrix, kernel, primal, dual = best_data
        certificate = verify_certificate(kernel, atoms, pairs, primal, dual)
        if Fraction(certificate["objective"]) != expected_minima[k]:
            raise AssertionError(certificate)
        orders.append(
            {
                "k": k,
                "antipodal_type_count": len(types),
                "ternary_atom_count": len(atoms),
                "full_sign_matrix_count": matrix_count,
                "objective_histogram": histogram,
                "minimum_L": expected_minima[k],
                "minimizing_R": matrix.tolist(),
                "minimizing_R_upper_triangle": upper_triangle(matrix),
                "kernel_diagonal_range": int(np.ptp(np.diag(kernel))),
                "certificate": certificate,
                "universal_range_bound": range_lower_bound(k),
            }
        )

    # A generous finite audit of the exact split bound. The report gives an
    # analytic all-k proof; this loop catches factor/rounding errors.
    range_audit = [range_lower_bound(k) | {"k": k} for k in range(2, 257)]
    if not all(row["strictly_exceeds_k_to_three_halves"] for row in range_audit):
        raise AssertionError("range lower bound failed")
    minimum_range_ratio = min(
        range_audit, key=lambda row: row["ratio_to_k_to_three_halves"]
    )
    range_audit_canonical = json.dumps(
        range_audit, sort_keys=True, separators=(",", ":")
    )

    payload = {
        "schema": "quadratic-signing-finite-fibre-cut-decomposition-v1",
        "classification": (
            "exact rational certificates for each saved LP optimum; exhaustive "
            "floating-point LP comparison over all R for k=2,3,4"
        ),
        "definition": (
            "L(R)=min sum_h |c_h| subject to K_R=sum_h c_h h h^T, "
            "h in {-1,0,1}^{2^(k-1)} modulo sign"
        ),
        "orders": orders,
        "range_bound_audit_through": 256,
        "range_bound_minimum_ratio": minimum_range_ratio,
        "range_bound_audit_sha256": hashlib.sha256(
            range_audit_canonical.encode()
        ).hexdigest(),
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    payload["canonical_payload_sha256"] = hashlib.sha256(canonical.encode()).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "output": str(args.output),
        "sha256": payload["canonical_payload_sha256"],
        "minimum_L": {str(row["k"]): row["minimum_L"] for row in orders},
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
