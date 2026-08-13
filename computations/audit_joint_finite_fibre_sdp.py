#!/usr/bin/env python3
"""Numerically audit the joint finite-fibre signed-elliptope certificate.

The proved certificate in ``artifacts/joint_finite_fibre_action_audit.md`` is

    Gamma(A) = min sum_i y_i
               subject to diag(y) +/- A/2 positive semidefinite.

This script verifies the saved signing witnesses and solves that small SDP for
orders 3 through 14.  The SDP values are numerical diagnostics; the universal
conference floor and the lift inequality are proved analytically in the
artifact and do not depend on the solver.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from pathlib import Path

import cvxpy as cp
import numpy as np


SOURCES: list[tuple[str, int, str, str]] = [
    *[
        (f"exact_m{n}", n, f"computations/results/exact_m{n}.json", "matrix")
        for n in range(3, 11)
    ],
    (
        "m8_orbit_class_1",
        8,
        "computations/results/m8_minimizer_orbits.json",
        "classes.1.representative_matrix",
    ),
    ("nested_m11", 11, "computations/results/nested_10_in_11_cap17.json", "matrix"),
    (
        "nested_m12",
        12,
        "computations/results/extension_nested_m11_to_12.json",
        "parent_matrix",
    ),
    (
        "bridge_m13",
        13,
        "computations/results/bridge_6_7_sign1_cap20.json",
        "parent_matrix",
    ),
    (
        "conference_m14",
        14,
        "computations/results/heuristic_m14_from_conference.json",
        "matrix",
    ),
]

KNOWN_VALUES = {
    3: 3,
    4: 4,
    5: 4,
    6: 5,
    7: 9,
    8: 10,
    9: 12,
    10: 13,
    11: 17,
    12: 18,
    13: 20,
    14: 21,
}


def load_nested(payload: dict[str, object], key: str) -> object:
    value: object = payload
    for part in key.split("."):
        if isinstance(value, dict):
            value = value[part]
        elif isinstance(value, list):
            value = value[int(part)]
        else:
            raise AssertionError((key, part))
    return value


def matrix_sha256(matrix: np.ndarray) -> str:
    canonical = json.dumps(matrix.astype(int).tolist(), separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest()


def exact_cap(matrix: np.ndarray) -> int:
    n = len(matrix)
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    signs = np.asarray([matrix[i, j] for i, j in edges], dtype=np.int64)
    best = 0
    for tail in itertools.product((-1, 1), repeat=n - 1):
        spin = (1,) + tail
        energy = sum(
            int(sign) * spin[i] * spin[j]
            for sign, (i, j) in zip(signs, edges)
        )
        best = max(best, abs(energy))
    return best


def solve_gamma(matrix: np.ndarray) -> dict[str, object]:
    n = len(matrix)
    y = cp.Variable(n)
    minus = cp.diag(y) - matrix / 2
    plus = cp.diag(y) + matrix / 2
    problem = cp.Problem(cp.Minimize(cp.sum(y)), [minus >> 0, plus >> 0])
    objective = problem.solve(
        solver="CLARABEL",
        tol_gap_abs=1e-8,
        tol_gap_rel=1e-8,
        tol_feas=1e-8,
        max_iter=1000,
    )
    if problem.status not in {"optimal", "optimal_inaccurate"}:
        raise AssertionError(problem.status)
    diagonal = np.asarray(y.value, dtype=np.float64)
    minimum_minus = float(np.min(np.linalg.eigvalsh(np.diag(diagonal) - matrix / 2)))
    minimum_plus = float(np.min(np.linalg.eigvalsh(np.diag(diagonal) + matrix / 2)))
    return {
        "solver": "CLARABEL",
        "status": problem.status,
        "objective": float(objective),
        "diagonal": diagonal.tolist(),
        "minimum_eigenvalue_diag_minus_A_over_2": minimum_minus,
        "minimum_eigenvalue_diag_plus_A_over_2": minimum_plus,
    }


def solve_correlation_endpoints(matrix: np.ndarray) -> dict[str, float | str]:
    """Ordinary degree-two endpoints, for normalization comparison only."""

    n = len(matrix)
    values = []
    statuses = []
    for sign in (1, -1):
        gram = cp.Variable((n, n), symmetric=True)
        problem = cp.Problem(
            cp.Maximize(sign * cp.trace(matrix @ gram) / 2),
            [gram >> 0, cp.diag(gram) == 1],
        )
        values.append(
            float(
                problem.solve(
                    solver="CLARABEL",
                    tol_gap_abs=1e-8,
                    tol_gap_rel=1e-8,
                    tol_feas=1e-8,
                    max_iter=1000,
                )
            )
        )
        statuses.append(problem.status)
    return {
        "positive_endpoint_U1": values[0],
        "negative_endpoint_minus_L1": values[1],
        "absolute_endpoint_R1": max(values),
        "centered_half_width_W1": sum(values) / 2,
        "statuses": ",".join(statuses),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    records = []
    for label, n, path_text, key in SOURCES:
        path = Path(path_text)
        payload = json.loads(path.read_text())
        matrix = np.asarray(load_nested(payload, key), dtype=np.float64)
        if matrix.shape != (n, n):
            raise AssertionError((n, matrix.shape))
        if not np.array_equal(matrix, matrix.T) or np.any(np.diag(matrix)):
            raise AssertionError((n, "not symmetric zero diagonal"))
        off_diagonal = matrix[~np.eye(n, dtype=bool)]
        if not np.all(np.isin(off_diagonal, (-1, 1))):
            raise AssertionError((n, "not a signing"))
        cap = exact_cap(matrix.astype(np.int64))
        if cap != KNOWN_VALUES[n]:
            raise AssertionError((n, cap, KNOWN_VALUES[n]))
        gamma = solve_gamma(matrix)
        correlation = solve_correlation_endpoints(matrix)
        floor = n * np.sqrt(n - 1) / 2
        if float(gamma["objective"]) < floor - 2e-7:
            raise AssertionError((n, gamma["objective"], floor))
        records.append(
            {
                "n": n,
                "label": label,
                "source": path_text,
                "source_key": key,
                "matrix_sha256": matrix_sha256(matrix),
                "certified_cap": cap,
                "cap_over_n_to_three_halves": cap / n**1.5,
                "conference_floor": floor,
                "gamma": gamma,
                "ordinary_correlation_sdp": correlation,
                "gamma_minus_cap": float(gamma["objective"]) - cap,
                "gamma_over_n_to_three_halves": float(gamma["objective"])
                / n**1.5,
                "exact_conference_identity": bool(
                    np.array_equal(
                        matrix.astype(np.int64) @ matrix.astype(np.int64),
                        (n - 1) * np.eye(n, dtype=np.int64),
                    )
                ),
            }
        )

    # The order-four Sylvester matrix is the smallest symmetric Hadamard
    # microkernel used in the theorem.
    hadamard4 = np.asarray(
        [
            [1, 1, 1, 1],
            [1, -1, 1, -1],
            [1, 1, -1, -1],
            [1, -1, -1, 1],
        ],
        dtype=np.int64,
    )
    if not np.array_equal(hadamard4 @ hadamard4, 4 * np.eye(4, dtype=np.int64)):
        raise AssertionError("invalid symmetric Hadamard kernel")

    output = {
        "schema": "quadratic-signing-joint-finite-fibre-sdp-audit-v1",
        "classification": (
            "exact cap and matrix checks; numerical SDP diagnostics only; "
            "analytic lift theorem and conference floor are in the artifact"
        ),
        "gamma_definition": (
            "min sum(y) subject to diag(y)-A/2 PSD and diag(y)+A/2 PSD"
        ),
        "symmetric_hadamard4": hadamard4.tolist(),
        "records": records,
    }
    canonical = json.dumps(output, sort_keys=True, separators=(",", ":"))
    output["canonical_payload_sha256"] = hashlib.sha256(canonical.encode()).hexdigest()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {
                "output": str(args.output),
                "canonical_payload_sha256": output["canonical_payload_sha256"],
                "cases": [row["label"] for row in records],
                "conference_orders": [
                    row["n"] for row in records if row["exact_conference_identity"]
                ],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
