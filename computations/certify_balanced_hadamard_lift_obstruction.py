#!/usr/bin/env python3
"""Certify a scalable obstruction for a balanced-diagonal Hadamard lift.

For the saved order-14 conference signing A, a balanced diagonal D, and the
symmetric Sylvester matrix H_k, define

    S(k) = A tensor H_k + D tensor (H_k - diag(H_k)).

The script verifies a Boolean witness at k=4 and its exact tensor
amplification for k=4^r.  No heuristic optimization is used by the
certificate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


DIAGONAL = np.asarray(
    [1, -1, -1, -1, 1, 1, -1, 1, 1, 1, -1, -1, 1, -1],
    dtype=np.int64,
)

BASE_WITNESS = np.asarray(
    [
        -1, -1, 1, 1, 1, -1, -1, -1, -1, 1, 1, 1, 1, 1,
        1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, 1, -1, 1,
        1, 1, 1, -1, 1, 1, 1, -1, 1, 1, 1, -1, -1, 1,
        1, 1, 1, 1, -1, -1, -1, -1, -1, 1, 1, 1, -1, -1,
    ],
    dtype=np.int64,
)

# Boolean +2 eigenvector of the order-four Sylvester matrix.
MICRO_EIGENVECTOR = np.asarray([-1, -1, -1, 1], dtype=np.int64)


def sylvester(order: int) -> np.ndarray:
    if order < 1 or order & (order - 1):
        raise ValueError("order must be a power of two")
    matrix = np.ones((1, 1), dtype=np.int64)
    while len(matrix) < order:
        matrix = np.block([[matrix, matrix], [matrix, -matrix]])
    return matrix


def lift(base: np.ndarray, order: int) -> np.ndarray:
    hadamard = sylvester(order)
    off_hadamard = hadamard - np.diag(np.diag(hadamard))
    return np.kron(base, hadamard) + np.kron(np.diag(DIAGONAL), off_hadamard)


def energy(matrix: np.ndarray, spin: np.ndarray) -> int:
    value = int(spin @ matrix @ spin)
    if value % 2:
        raise AssertionError("quadratic value must be even")
    return value // 2


def matrix_hash(matrix: np.ndarray) -> str:
    return hashlib.sha256(matrix.astype(np.int8).tobytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        type=Path,
        default=Path("computations/results/heuristic_m14_from_conference.json"),
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    source = json.loads(args.source.read_text())
    base = np.asarray(source["matrix"], dtype=np.int64)
    n = len(base)
    if base.shape != (14, 14):
        raise AssertionError("expected the saved order-14 signing")
    if not np.array_equal(base, base.T) or np.any(np.diag(base)):
        raise AssertionError("base is not a symmetric zero-diagonal signing")
    if not np.array_equal(base @ base, 13 * np.eye(14, dtype=np.int64)):
        raise AssertionError("base conference identity failed")
    if int(source["profile"]["M"]) != 21:
        raise AssertionError("saved exact cap changed")
    if int(DIAGONAL.sum()) != 0:
        raise AssertionError("macro diagonal must be balanced")

    h4 = sylvester(4)
    if not np.array_equal(h4 @ MICRO_EIGENVECTOR, 2 * MICRO_EIGENVECTOR):
        raise AssertionError("micro Boolean eigenvector failed")
    s4 = lift(base, 4)
    if s4.shape != (56, 56) or np.any(np.diag(s4)):
        raise AssertionError("invalid order-56 lift")
    off = s4[~np.eye(56, dtype=bool)]
    if not np.all(np.abs(off) == 1):
        raise AssertionError("lift has a non-sign off-diagonal entry")
    base_energy = energy(s4, BASE_WITNESS)
    if base_energy != 220:
        raise AssertionError(f"base witness changed: {base_energy}")

    records = []
    tensor_witness = BASE_WITNESS.copy()
    for exponent in range(1, 5):
        micro_order = 4**exponent
        if exponent > 1:
            tensor_witness = np.kron(tensor_witness, MICRO_EIGENVECTOR)
        total_order = 14 * micro_order
        expected = 220 * (4 ** (exponent - 1)) ** 1.5
        expected_integer = 220 * 8 ** (exponent - 1)
        if int(expected) != expected_integer:
            raise AssertionError("scaling arithmetic failed")
        # Materialize modest cases as an independent exact audit.  The general
        # identity follows from Kronecker multiplication.
        checked_energy = None
        if exponent <= 3:
            matrix = lift(base, micro_order)
            checked_energy = energy(matrix, tensor_witness)
            if checked_energy != expected_integer:
                raise AssertionError((exponent, checked_energy, expected_integer))
        records.append(
            {
                "exponent": exponent,
                "micro_order": micro_order,
                "total_order": total_order,
                "certified_energy_lower_bound": expected_integer,
                "materialized_energy_check": checked_energy,
                "normalized_lower_bound": expected_integer / total_order**1.5,
            }
        )

    normalized = 220 / 56**1.5
    payload = {
        "schema": "quadratic-signing-balanced-hadamard-lift-obstruction-v1",
        "classification": (
            "exact fixed witness, exact Kronecker amplification, and proved "
            "scalable cap lower bound; no heuristic cap claim"
        ),
        "source": str(args.source),
        "source_matrix_sha256": matrix_hash(base),
        "macro_diagonal": DIAGONAL.tolist(),
        "macro_diagonal_trace": int(DIAGONAL.sum()),
        "definition": (
            "S_r=A tensor H_(4^r) + D tensor "
            "(H_(4^r)-diag(H_(4^r)))"
        ),
        "order_56_witness": BASE_WITNESS.tolist(),
        "order_56_witness_energy": base_energy,
        "order_56_matrix_sha256": matrix_hash(s4),
        "micro_boolean_eigenvector": MICRO_EIGENVECTOR.tolist(),
        "normalized_cap_lower_bound": normalized,
        "conference_constant": 0.5,
        "strict_constant_gap": normalized - 0.5,
        "records": records,
    }
    print(
        "verified balanced Hadamard lift obstruction: "
        f"cap(S_r)/N_r^(3/2) >= {normalized:.15f} > 1/2"
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
