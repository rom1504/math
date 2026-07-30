#!/usr/bin/env python3
"""Find and verify a Boolean Hadamard eigenvector for a doubled conference signing.

If H=A+diag(-I_q,I_q) is the symmetric Hadamard completion and its order
N=2q is a square, a sign vector satisfying Hx=+/-sqrt(N)x attains the spectral
upper bound.  The saved vector and matrix multiplication are an independently
checkable proof that the signing cap is N^(3/2)/2.  CP-SAT is used only to
discover the vector.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from ortools.sat.python import cp_model

from exact_mn_milp import stable_matrix_hash


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("construction", type=Path)
    parser.add_argument("--time-limit", type=float, default=300.0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = json.loads(args.construction.read_text())
    a = np.asarray(payload["parent_matrix"], dtype=np.int64)
    n = len(a)
    q = n // 2
    root = math.isqrt(n)
    if root * root != n:
        raise ValueError("parent order is not a square")
    correction = np.diag([-1] * q + [1] * q).astype(np.int64)
    h = a + correction
    if not np.array_equal(h @ h.T, n * np.eye(n, dtype=np.int64)):
        raise AssertionError("Hadamard completion failed")

    certificates = []
    for eigenvalue in (root, -root):
        model = cp_model.CpModel()
        bits = [model.new_bool_var(f"z_{i}") for i in range(n)]
        model.add(bits[0] == 0)
        row_sums = h.sum(axis=1, dtype=np.int64)
        for i in range(n):
            lhs = int(row_sums[i]) - 2 * sum(
                int(h[i, j]) * bits[j] for j in range(n)
            )
            model.add(lhs == eigenvalue * (1 - 2 * bits[i]))
        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = args.time_limit
        solver.parameters.num_search_workers = 1
        status = solver.solve(model)
        status_name = solver.status_name(status)
        record: dict[str, object] = {
            "eigenvalue": eigenvalue,
            "solver_status": status_name,
            "solver_wall_time_seconds": solver.wall_time,
        }
        if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            x = np.asarray([1 - 2 * solver.value(bit) for bit in bits], dtype=np.int64)
            if not np.array_equal(h @ x, eigenvalue * x):
                raise AssertionError("saved vector is not the claimed eigenvector")
            energy = int(x @ a @ x // 2)
            if energy != eigenvalue * n // 2:
                raise AssertionError((energy, eigenvalue * n // 2))
            record["spins"] = [int(value) for value in x]
            record["verified_energy"] = energy
        certificates.append(record)
        print(
            f"eigenvalue={eigenvalue:+d} status={status_name} "
            f"energy={record.get('verified_energy')} wall={solver.wall_time:.6f}s"
        )
    witnessed = any("spins" in record for record in certificates)
    output = {
        "schema": "quadratic-signing-conference-double-eigen-v1",
        "classification": (
            "proved exact fixed-signing cap from an explicit Boolean eigenvector "
            "and the Hadamard spectral bound"
            if witnessed
            else "inconclusive Boolean-eigenvector search; no cap improvement"
        ),
        "source": str(args.construction),
        "n": n,
        "matrix_sha256": stable_matrix_hash(a),
        "hadamard_eigenvalue_magnitude": root,
        "certified_cap": n * root // 2 if witnessed else None,
        "normalized_cap": 0.5 if witnessed else None,
        "certificates": certificates,
        "proof": (
            [
                "H=A+diag(-I_q,I_q) and H H^T=nI",
                "x^T diag(-I_q,I_q) x=0 for every sign vector x",
                "therefore |x^T A x|/2<=n^(3/2)/2",
                "the saved Boolean eigenvector attains equality",
            ]
            if witnessed
            else []
        ),
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
        print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
