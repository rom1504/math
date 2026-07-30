#!/usr/bin/env python3
"""Exhaustively optimize a one-vertex extension of a fixed signing.

For a fixed n-vertex matrix A and a sign row r, maximizing over the new spin
gives the exact identity

    M([[A,r],[r^T,0]]) = max_x (|H_A(x)| + |r dot x|).

The choices r and -r are equivalent, so r[0]=+1 is fixed.  This script
enumerates the remaining 2^(n-1) rows, saves every best row (subject to an
optional output limit), and exhaustively verifies the assembled witness.
The result is a proved finite statement about extensions of the supplied
matrix, not about all order-(n+1) signings.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from exact_mn_milp import exact_profile, projective_spins, stable_matrix_hash


def load_matrix(path: Path) -> np.ndarray:
    payload = json.loads(path.read_text())
    key = "matrix" if "matrix" in payload else "parent_matrix"
    matrix = np.asarray(payload[key], dtype=np.int8)
    if not np.array_equal(matrix, matrix.T) or np.any(np.diag(matrix)):
        raise ValueError(f"invalid signing matrix in {path}")
    return matrix


def assemble(matrix: np.ndarray, row: np.ndarray) -> np.ndarray:
    n = len(matrix)
    parent = np.zeros((n + 1, n + 1), dtype=np.int8)
    parent[:n, :n] = matrix
    parent[:n, n] = row
    parent[n, :n] = row
    return parent


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("child", type=Path)
    parser.add_argument("--row-output-limit", type=int, default=100)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    matrix = load_matrix(args.child)
    n = len(matrix)
    spins = projective_spins(n).astype(np.int64)
    energies = np.einsum("bi,ij,bj->b", spins, matrix.astype(np.int64), spins) // 2
    candidate_rows = projective_spins(n).astype(np.int64)
    best_cap: int | None = None
    best_rows: list[np.ndarray] = []
    histogram: dict[int, int] = {}
    for row in candidate_rows:
        cap = int(np.max(np.abs(energies) + np.abs(spins @ row)))
        histogram[cap] = histogram.get(cap, 0) + 1
        if best_cap is None or cap < best_cap:
            best_cap = cap
            best_rows = [row.copy()]
        elif cap == best_cap and len(best_rows) < args.row_output_limit:
            best_rows.append(row.copy())

    assert best_cap is not None and best_rows
    parent = assemble(matrix, best_rows[0])
    profile = exact_profile(parent)
    if profile["M"] != best_cap:
        raise AssertionError((profile["M"], best_cap))
    payload = {
        "schema": "quadratic-signing-exact-vertex-extension-v1",
        "classification": "proved exhaustive finite optimization for extensions of the supplied child",
        "child": str(args.child),
        "child_order": n,
        "candidate_rows_checked": len(candidate_rows),
        "best_extension_cap": best_cap,
        "best_row_count_saved": len(best_rows),
        "best_rows": [[int(v) for v in row] for row in best_rows],
        "extension_cap_histogram": {
            str(cap): count for cap, count in sorted(histogram.items())
        },
        "parent_matrix": [[int(v) for v in row] for row in parent],
        "parent_matrix_sha256": stable_matrix_hash(parent),
        "parent_profile": profile,
    }
    print(
        f"child n={n}: checked={len(candidate_rows)} best_extension_cap={best_cap} "
        f"saved_rows={len(best_rows)} hash={payload['parent_matrix_sha256']}"
    )
    print(f"extension cap histogram: {payload['extension_cap_histogram']}")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
