#!/usr/bin/env python3
"""Construct the scalable doubled signing from a Paley conference matrix.

For a prime p=1 mod 4, the symmetric Paley conference matrix S has order
q=p+1 and S^2=pI.  The order-2q signing

    A = [[S, S+I], [S+I, -S]]

has zero diagonal and sign off-diagonal entries.  It satisfies

    (A^2-(2q-1)I)^2 = 4(q-1)I.

The matrix identities are exact and scalable.  This script does not claim a
useful Boolean cap; that is measured separately.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from exact_mn_milp import exact_profile, stable_matrix_hash


def is_prime(value: int) -> bool:
    if value < 2:
        return False
    divisor = 2
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 1
    return True


def legendre(value: int, prime: int) -> int:
    value %= prime
    if value == 0:
        return 0
    residue = pow(value, (prime - 1) // 2, prime)
    return 1 if residue == 1 else -1


def paley_conference(prime: int) -> np.ndarray:
    if not is_prime(prime) or prime % 4 != 1:
        raise ValueError("prime must be a prime congruent to 1 modulo 4")
    order = prime + 1
    matrix = np.zeros((order, order), dtype=np.int8)
    matrix[0, 1:] = matrix[1:, 0] = 1
    for i in range(prime):
        for j in range(prime):
            if i != j:
                matrix[i + 1, j + 1] = legendre(i - j, prime)
    target = prime * np.eye(order, dtype=np.int64)
    if not np.array_equal(matrix.astype(np.int64) @ matrix.astype(np.int64), target):
        raise AssertionError("Paley conference identity failed")
    return matrix


def double_conference(matrix: np.ndarray) -> np.ndarray:
    order = len(matrix)
    identity = np.eye(order, dtype=np.int8)
    return np.block([[matrix, matrix + identity], [matrix + identity, -matrix]])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("prime", type=int)
    parser.add_argument("--exhaustive-limit", type=int, default=22)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    conference = paley_conference(args.prime)
    parent = double_conference(conference)
    q = len(conference)
    identity = np.eye(2 * q, dtype=np.int64)
    parent64 = parent.astype(np.int64)
    quartic = (parent64 @ parent64 - (2 * q - 1) * identity)
    if not np.array_equal(quartic @ quartic, 4 * (q - 1) * identity):
        raise AssertionError("doubled quartic identity failed")
    diagonal_correction = np.diag([-1] * q + [1] * q).astype(np.int64)
    hadamard = parent64 + diagonal_correction
    if not np.array_equal(hadamard @ hadamard.T, 2 * q * identity):
        raise AssertionError("symmetric Hadamard identity failed")
    payload: dict[str, object] = {
        "schema": "quadratic-signing-conference-double-v1",
        "classification": "proved explicit construction and exact integer matrix identities; Boolean cap open unless profile is present",
        "paley_prime": args.prime,
        "conference_order": q,
        "parent_order": 2 * q,
        "conference_matrix": [[int(v) for v in row] for row in conference],
        "conference_matrix_sha256": stable_matrix_hash(conference),
        "parent_matrix": [[int(v) for v in row] for row in parent],
        "parent_matrix_sha256": stable_matrix_hash(parent),
        "symmetric_hadamard_matrix": [
            [int(v) for v in row] for row in hadamard
        ],
        "identities": {
            "conference_square": f"S^2={q-1}I",
            "parent_quartic": f"(A^2-{2*q-1}I)^2={4*(q-1)}I",
            "hadamard_completion": (
                "H=A+diag(-I_q,I_q), H H^T=2q I, and "
                "x^T A x=x^T H x for every sign vector x"
            ),
            "verified_exact_integer_arithmetic": True,
        },
        "exact_boolean_identity": (
            "for x in signs and J subset [q], parent energy equals "
            "2 H_S(x)-4 H_(S[J])(x_J)+q-2|J|"
        ),
    }
    if q <= args.exhaustive_limit:
        payload["conference_profile"] = exact_profile(conference)
    if 2 * q <= args.exhaustive_limit:
        payload["parent_profile"] = exact_profile(parent)
        payload["classification"] = (
            "proved explicit construction, exact integer identities, and exhaustive Boolean profile"
        )
    elif q <= args.exhaustive_limit:
        payload["classification"] = (
            "proved explicit construction and exact integer identities, with an "
            "exhaustive Boolean profile for the conference child only"
        )
    print(
        f"Paley p={args.prime}: conference_order={q} parent_order={2*q} "
        f"parent_hash={payload['parent_matrix_sha256']}"
    )
    if "parent_profile" in payload:
        print(f"exhaustive parent cap={payload['parent_profile']['M']}")
    elif "conference_profile" in payload:
        print(f"exhaustive conference cap={payload['conference_profile']['M']}")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
