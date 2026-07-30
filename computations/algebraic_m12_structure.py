#!/usr/bin/env python3
"""Verify the exact 6+6 algebraic structure of the saved M_12 witness.

The script discovers the two support components of

    D = (A^2 - 11 I) / 2

and verifies, after the induced permutation, that

    A = [[S, C], [C^T, T]],
    S^2 = T^2 = 5 I,
    S C + C T = 0,
    C C^T = 6 I + 2 S,
    C^T C = 6 I - 2 T.

Consequently (A^2-11I)^2=20I.  Both diagonal blocks are exact order-6
minimizers, and the cap of the assembled matrix is recomputed exhaustively.
All claims in the output are exact integer identities about this witness.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from exact_mn_milp import exact_profile, stable_matrix_hash


def load_matrix(path: Path) -> np.ndarray:
    payload = json.loads(path.read_text())
    key = "matrix" if "matrix" in payload else "parent_matrix"
    matrix = np.asarray(payload[key], dtype=np.int64)
    if matrix.shape != (12, 12):
        raise ValueError("this certificate expects an order-12 witness")
    return matrix


def support_components(matrix: np.ndarray) -> list[list[int]]:
    support = np.abs(matrix)
    seen: set[int] = set()
    components: list[list[int]] = []
    for start in range(len(matrix)):
        if start in seen:
            continue
        seen.add(start)
        stack = [start]
        component = []
        while stack:
            vertex = stack.pop()
            component.append(vertex)
            for neighbor in np.flatnonzero(support[vertex]):
                item = int(neighbor)
                if item not in seen:
                    seen.add(item)
                    stack.append(item)
        components.append(sorted(component))
    return components


def zero(matrix: np.ndarray) -> bool:
    return bool(np.all(matrix == 0))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "witness",
        type=Path,
        default=Path("computations/results/extension_nested_m11_to_12.json"),
        nargs="?",
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    a = load_matrix(args.witness)
    identity_12 = np.eye(12, dtype=np.int64)
    remainder = a @ a - 11 * identity_12
    if np.any(remainder % 2):
        raise AssertionError("A^2-11I is not even")
    d = remainder // 2
    components = support_components(d)
    if sorted(map(len, components)) != [6, 6]:
        raise AssertionError(components)
    left, right = components
    s = a[np.ix_(left, left)]
    t = a[np.ix_(right, right)]
    c = a[np.ix_(left, right)]
    identity_6 = np.eye(6, dtype=np.int64)
    checks = {
        "S_squared_equals_5I": zero(s @ s - 5 * identity_6),
        "T_squared_equals_5I": zero(t @ t - 5 * identity_6),
        "intertwining_SC_plus_CT_zero": zero(s @ c + c @ t),
        "CCt_equals_6I_plus_2S": zero(c @ c.T - 6 * identity_6 - 2 * s),
        "CtC_equals_6I_minus_2T": zero(c.T @ c - 6 * identity_6 + 2 * t),
        "quartic_parent_identity": zero(
            np.linalg.matrix_power(a, 4) - 22 * (a @ a) + 101 * identity_12
        ),
    }
    if not all(checks.values()):
        raise AssertionError(checks)
    parent_profile = exact_profile(a)
    left_profile = exact_profile(s)
    right_profile = exact_profile(t)
    if (parent_profile["M"], left_profile["M"], right_profile["M"]) != (18, 5, 5):
        raise AssertionError(
            (parent_profile["M"], left_profile["M"], right_profile["M"])
        )
    singular_values = np.linalg.svd(c.astype(float), compute_uv=False)
    payload = {
        "schema": "quadratic-signing-algebraic-m12-structure-v1",
        "classification": "proved exact integer identities and exhaustive finite energy profiles for the supplied witness",
        "source": str(args.witness),
        "matrix_sha256": stable_matrix_hash(a),
        "partition": [left, right],
        "checks": checks,
        "profiles": {
            "parent": parent_profile,
            "left_child": left_profile,
            "right_child": right_profile,
        },
        "bridge": [[int(value) for value in row] for row in c],
        "bridge_singular_values": [float(value) for value in singular_values],
        "deductions": [
            "D=(A^2-11I)/2 is block diagonal on two six-vertex components",
            "the diagonal blocks S and T are symmetric conference matrices",
            "the parent characteristic polynomial is (lambda^4-22 lambda^2+101)^3",
            "the bridge singular values are sqrt(6+2sqrt(5)) and sqrt(6-2sqrt(5)), each with multiplicity three",
        ],
    }
    print("verified exact 6+6 conference/bridge structure of the M_12 witness")
    print(f"partition={components} parent_M=18 child_caps=5,5")
    print(f"checks={checks}")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
