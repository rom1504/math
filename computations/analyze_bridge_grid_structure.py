#!/usr/bin/env python3
"""Measure algebraic structure in the feasible fixed-child bridge grid.

For each saved feasible bridge, this computes exact intertwining and Gram
residuals, spectra, and parent fourth moments.  The affine Gram fit asks how
close CC^T is to alpha*I+beta*A (and similarly on the right).  Results concern
the saved solver witnesses only; CP-SAT does not optimize these secondary
statistics.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def affine_fit(target: np.ndarray, child: np.ndarray) -> dict[str, object]:
    n = len(child)
    identity = np.eye(n, dtype=float)
    design = np.column_stack((identity.ravel(), child.astype(float).ravel()))
    coefficients, *_ = np.linalg.lstsq(design, target.astype(float).ravel(), rcond=None)
    fitted = coefficients[0] * identity + coefficients[1] * child
    residual = target - fitted
    return {
        "alpha": float(coefficients[0]),
        "beta": float(coefficients[1]),
        "residual_frobenius": float(np.linalg.norm(residual)),
        "target_frobenius": float(np.linalg.norm(target)),
        "relative_residual": float(np.linalg.norm(residual) / np.linalg.norm(target)),
        "exact_affine": bool(np.max(np.abs(residual)) < 1e-9),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "grid",
        type=Path,
        default=Path("computations/results/bridge_grid_through_12.json"),
        nargs="?",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    grid = json.loads(args.grid.read_text())
    records = []
    for row in grid["rows"]:
        if row["status"] != "OPTIMAL":
            continue
        payload = json.loads(Path(row["result"]).read_text())
        a_payload = json.loads(Path(payload["child_a"]).read_text())
        b_payload = json.loads(Path(payload["child_b"]).read_text())
        a = np.asarray(a_payload["matrix"], dtype=np.int64)
        b = payload["sign_b"] * np.asarray(b_payload["matrix"], dtype=np.int64)
        c = np.asarray(payload["bridge"], dtype=np.int64)
        parent = np.asarray(payload["parent_matrix"], dtype=np.int64)
        intertwining = a @ c + c @ b
        square = parent @ parent
        eigenvalues = np.linalg.eigvalsh(parent.astype(float))
        distinct = []
        for value in eigenvalues:
            if not distinct or abs(value - distinct[-1]) > 1e-7:
                distinct.append(float(value))
        records.append(
            {
                "m": row["m"],
                "n": row["n"],
                "sign_b": row["sign_b"],
                "parent_M": payload["parent_profile"]["M"],
                "parent_operator_norm": float(np.max(np.abs(eigenvalues))),
                "parent_trace_A4": int(np.trace(square @ square)),
                "parent_distinct_eigenvalue_count": len(distinct),
                "parent_distinct_eigenvalues": distinct,
                "bridge_singular_values": [
                    float(v) for v in np.linalg.svd(c.astype(float), compute_uv=False)
                ],
                "intertwining_frobenius": float(np.linalg.norm(intertwining)),
                "intertwining_max_absolute": int(np.max(np.abs(intertwining))),
                "left_gram_affine_fit": affine_fit(c @ c.T, a),
                "right_gram_affine_fit": affine_fit(c.T @ c, b),
                "source": row["result"],
            }
        )
    output = {
        "schema": "quadratic-signing-bridge-grid-structure-v1",
        "classification": "exact statistics for saved feasible bridge witnesses; not optimized secondary invariants",
        "source": str(args.grid),
        "records": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    exact = [
        record
        for record in records
        if record["intertwining_max_absolute"] == 0
        and record["left_gram_affine_fit"]["exact_affine"]
        and record["right_gram_affine_fit"]["exact_affine"]
    ]
    print(f"analyzed {len(records)} feasible bridges; exact algebraic triples={len(exact)}")
    for record in exact:
        print(
            f"  {record['m']}+{record['n']} sign={record['sign_b']:+d} "
            f"eigenvalue_count={record['parent_distinct_eigenvalue_count']}"
        )
    ranked = sorted(records, key=lambda record: record["intertwining_frobenius"])
    print("smallest intertwining residuals:")
    for record in ranked[:10]:
        print(
            f"  {record['m']}+{record['n']} sign={record['sign_b']:+d} "
            f"residual={record['intertwining_frobenius']:.6f}"
        )
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
