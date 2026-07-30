#!/usr/bin/env python3
"""Exact MILP for the quadratic-signing minimum M_n.

The model fixes vertex 0 by Seidel switching, so a[0,j]=+1.  It uses one
binary variable z[i,j] for every remaining edge, with a[i,j]=1-2*z[i,j],
and an integer cap t.  For every projective Boolean spin (x[0]=+1), it adds

    -t <= sum_{i<j} a[i,j] x[i] x[j] <= t.

The optional basic symmetry constraints are valid under permutations of the
non-root vertices and global negation followed by re-gauging:

* vertex 1 has minimum negative degree in the rooted internal graph;
* its incident signs are sorted, positive before negative; and
* its first internal edge is positive.

SciPy delegates the branch-and-bound computation to HiGHS.  A zero-gap
optimal result is a solver-certified computation, not a formal proof object.
The saved JSON includes an independently recomputed exhaustive energy profile
and the solver's primal/dual bounds.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np
import scipy
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import coo_matrix, csr_matrix, vstack


@dataclass(frozen=True)
class ModelData:
    n: int
    internal_edges: tuple[tuple[int, int], ...]
    spins: np.ndarray
    pair_products: np.ndarray
    constraint: LinearConstraint
    bounds: Bounds
    integrality: np.ndarray
    objective: np.ndarray


def projective_spins(n: int) -> np.ndarray:
    """All Boolean spins with x[0]=+1, in deterministic bit order."""
    count = 1 << (n - 1)
    bits = np.arange(count, dtype=np.uint64)[:, None]
    shifts = np.arange(n - 1, dtype=np.uint64)[None, :]
    tail = 1 - 2 * ((bits >> shifts) & 1).astype(np.int8)
    return np.column_stack((np.ones(count, dtype=np.int8), tail))


def build_model(n: int, symmetry: str = "basic") -> ModelData:
    if n < 3:
        raise ValueError("the implemented rooted symmetry model requires n>=3")
    if symmetry not in {"none", "basic"}:
        raise ValueError(f"unknown symmetry mode: {symmetry}")

    edges = tuple(combinations(range(1, n), 2))
    edge_index = {edge: k for k, edge in enumerate(edges)}
    spins = projective_spins(n)
    pair_products = np.column_stack(
        [spins[:, i] * spins[:, j] for i, j in edges]
    ).astype(np.int8, copy=False)

    root_part = spins[:, 1:].sum(axis=1, dtype=np.int32)
    internal_all_positive = pair_products.sum(axis=1, dtype=np.int32)
    constant = root_part + internal_all_positive
    coeff = (-2 * pair_products).astype(np.float64)
    cap_column = -np.ones((spins.shape[0], 1), dtype=np.float64)
    energy_upper = csr_matrix(np.hstack((coeff, cap_column)))
    energy_lower = csr_matrix(np.hstack((-coeff, cap_column)))
    matrix_parts: list[csr_matrix] = [energy_upper, energy_lower]
    lower_parts = [np.full(spins.shape[0], -np.inf)] * 2
    upper_parts = [-constant.astype(float), constant.astype(float)]

    variable_count = len(edges) + 1
    lower_bounds = np.zeros(variable_count, dtype=float)
    upper_bounds = np.ones(variable_count, dtype=float)
    upper_bounds[-1] = np.inf

    if symmetry == "basic":
        # Choose vertex 1 with minimum negative degree.  Permute its positive
        # neighbors before its negative neighbors.  Complement the rooted
        # internal graph if necessary so that the first edge is positive.
        upper_bounds[edge_index[(1, 2)]] = 0.0
        row_indices: list[int] = []
        col_indices: list[int] = []
        values: list[float] = []
        rhs: list[float] = []

        def add_row(terms: dict[int, float], bound: float = 0.0) -> None:
            row = len(rhs)
            for col, value in terms.items():
                if value:
                    row_indices.append(row)
                    col_indices.append(col)
                    values.append(value)
            rhs.append(bound)

        # z[1,j] is nondecreasing in j.
        for j in range(2, n - 1):
            add_row({edge_index[(1, j)]: 1.0, edge_index[(1, j + 1)]: -1.0})

        # deg_-(1) <= deg_-(i) for each i>1.
        incident_1 = [edge_index[(1, j)] for j in range(2, n)]
        for i in range(2, n):
            terms: dict[int, float] = {}
            for col in incident_1:
                terms[col] = terms.get(col, 0.0) + 1.0
            for j in range(1, n):
                if i == j:
                    continue
                edge = (j, i) if j < i else (i, j)
                col = edge_index[edge]
                terms[col] = terms.get(col, 0.0) - 1.0
            add_row(terms)

        if rhs:
            sym = coo_matrix(
                (values, (row_indices, col_indices)),
                shape=(len(rhs), variable_count),
            ).tocsr()
            matrix_parts.append(sym)
            lower_parts.append(np.full(len(rhs), -np.inf))
            upper_parts.append(np.asarray(rhs, dtype=float))

    matrix = vstack(matrix_parts, format="csr")
    lower = np.concatenate(lower_parts)
    upper = np.concatenate(upper_parts)
    constraint = LinearConstraint(matrix, lower, upper)
    bounds = Bounds(lower_bounds, upper_bounds)
    integrality = np.ones(variable_count, dtype=np.uint8)
    objective = np.zeros(variable_count, dtype=float)
    objective[-1] = 1.0
    return ModelData(
        n=n,
        internal_edges=edges,
        spins=spins,
        pair_products=pair_products,
        constraint=constraint,
        bounds=bounds,
        integrality=integrality,
        objective=objective,
    )


def matrix_from_solution(model: ModelData, solution: np.ndarray) -> np.ndarray:
    n = model.n
    matrix = np.zeros((n, n), dtype=np.int8)
    matrix[0, 1:] = 1
    matrix[1:, 0] = 1
    rounded = np.rint(solution[: len(model.internal_edges)]).astype(np.int8)
    for (i, j), z in zip(model.internal_edges, rounded):
        sign = 1 - 2 * int(z)
        matrix[i, j] = matrix[j, i] = sign
    return matrix


def exact_profile(matrix: np.ndarray) -> dict[str, object]:
    spins = projective_spins(matrix.shape[0]).astype(np.int64)
    energies = np.einsum("bi,ij,bj->b", spins, matrix.astype(np.int64), spins) // 2
    max_energy = int(energies.max())
    min_energy = int(energies.min())
    cap = int(np.abs(energies).max())
    top = spins[energies == max_energy]
    bottom = spins[energies == min_energy]
    row_sums = matrix.sum(axis=1, dtype=np.int64)
    eigenvalues = np.linalg.eigvalsh(matrix.astype(float))
    histogram_values, histogram_counts = np.unique(energies, return_counts=True)
    return {
        "M": cap,
        "P": max_energy,
        "Q": -min_energy,
        "projective_top_count": int(len(top)),
        "projective_bottom_count": int(len(bottom)),
        "row_sums": [int(v) for v in row_sums],
        "eigenvalues": [float(v) for v in eigenvalues],
        "energy_histogram": {
            str(int(v)): int(c) for v, c in zip(histogram_values, histogram_counts)
        },
        "top_spins": [[int(v) for v in row] for row in top],
        "bottom_spins": [[int(v) for v in row] for row in bottom],
    }


def stable_matrix_hash(matrix: np.ndarray) -> str:
    return hashlib.sha256(matrix.astype(np.int8).tobytes(order="C")).hexdigest()


def solve(args: argparse.Namespace) -> int:
    started = time.time()
    model = build_model(args.n, args.symmetry)
    print(
        f"model n={args.n} internal_binary={len(model.internal_edges)} "
        f"projective_spins={len(model.spins)} "
        f"linear_rows={model.constraint.A.shape[0]} symmetry={args.symmetry}",
        flush=True,
    )
    objective = model.objective.copy()
    lower_bounds = np.asarray(model.bounds.lb, dtype=float).copy()
    upper_bounds = np.asarray(model.bounds.ub, dtype=float).copy()
    if args.lower_bound is not None:
        lower_bounds[-1] = max(lower_bounds[-1], float(args.lower_bound))
    if args.upper_bound is not None:
        upper_bounds[-1] = min(upper_bounds[-1], float(args.upper_bound))
    if args.decision_cap is not None:
        lower_bounds[-1] = upper_bounds[-1] = float(args.decision_cap)
        objective[:] = 0.0
    if lower_bounds[-1] > upper_bounds[-1]:
        raise ValueError("inconsistent cap bounds")
    bounds = Bounds(lower_bounds, upper_bounds)

    options: dict[str, object] = {
        "disp": True,
        "presolve": True,
        "mip_rel_gap": 0.0,
    }
    if args.time_limit is not None:
        options["time_limit"] = float(args.time_limit)
    if args.node_limit is not None:
        options["node_limit"] = int(args.node_limit)

    result = milp(
        objective,
        integrality=model.integrality,
        bounds=bounds,
        constraints=model.constraint,
        options=options,
    )
    elapsed = time.time() - started
    print(f"status={result.status} success={result.success} message={result.message}")
    print(
        f"fun={result.fun} dual={getattr(result, 'mip_dual_bound', None)} "
        f"gap={getattr(result, 'mip_gap', None)} "
        f"nodes={getattr(result, 'mip_node_count', None)} elapsed={elapsed:.6f}s",
        flush=True,
    )

    payload: dict[str, object] = {
        "schema": "quadratic-signing-exact-milp-v1",
        "n": args.n,
        "normalization": "M_n=max_x |sum_{i<j} a_ij x_i x_j|",
        "symmetry": args.symmetry,
        "mode": "decision" if args.decision_cap is not None else "optimization",
        "cap_bounds": {
            "lower": None if args.lower_bound is None else args.lower_bound,
            "upper": None if args.upper_bound is None else args.upper_bound,
            "decision": args.decision_cap,
        },
        "model": {
            "internal_binary_variables": len(model.internal_edges),
            "total_integer_variables": len(model.objective),
            "projective_spin_constraints": 2 * len(model.spins),
            "total_linear_rows": int(model.constraint.A.shape[0]),
        },
        "solver": {
            "scipy_version": scipy.__version__,
            "python_version": platform.python_version(),
            "status": int(result.status),
            "success": bool(result.success),
            "message": str(result.message),
            "objective": None if result.fun is None else float(result.fun),
            "dual_bound": (
                None
                if getattr(result, "mip_dual_bound", None) is None
                else float(result.mip_dual_bound)
            ),
            "relative_gap": (
                None
                if getattr(result, "mip_gap", None) is None
                else float(result.mip_gap)
            ),
            "node_count": (
                None
                if getattr(result, "mip_node_count", None) is None
                else int(result.mip_node_count)
            ),
            "elapsed_seconds": elapsed,
            "time_limit_seconds": args.time_limit,
            "node_limit": args.node_limit,
        },
    }

    exit_code = 0 if result.x is not None else 2
    if result.x is not None:
        matrix = matrix_from_solution(model, result.x)
        profile = exact_profile(matrix)
        cap_variable = int(round(float(result.x[-1])))
        if profile["M"] > cap_variable:
            raise AssertionError((profile["M"], cap_variable))
        payload["matrix"] = [[int(v) for v in row] for row in matrix]
        payload["matrix_sha256"] = stable_matrix_hash(matrix)
        payload["profile"] = profile
        print(
            f"verified profile M={profile['M']} P={profile['P']} Q={profile['Q']} "
            f"sha256={payload['matrix_sha256']}",
            flush=True,
        )

    if args.output:
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(f"wrote {output}", flush=True)
    return exit_code


def verify(path: Path) -> int:
    payload = json.loads(path.read_text())
    matrix = np.asarray(payload["matrix"], dtype=np.int8)
    if stable_matrix_hash(matrix) != payload["matrix_sha256"]:
        raise AssertionError("matrix hash mismatch")
    if not np.array_equal(matrix, matrix.T):
        raise AssertionError("matrix is not symmetric")
    if np.any(np.diag(matrix)):
        raise AssertionError("matrix diagonal is not zero")
    if set(matrix[~np.eye(len(matrix), dtype=bool)]) != {-1, 1}:
        raise AssertionError("off-diagonal entries are not signs")
    profile = exact_profile(matrix)
    if profile != payload["profile"]:
        raise AssertionError("saved exhaustive profile does not recompute")
    print(
        f"verified {path}: n={len(matrix)} M={profile['M']} "
        f"P={profile['P']} Q={profile['Q']} hash={payload['matrix_sha256']}"
    )
    return 0


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    solve_parser = sub.add_parser("solve", help="build and solve one exact MILP")
    solve_parser.add_argument("n", type=int)
    solve_parser.add_argument("--symmetry", choices=("none", "basic"), default="basic")
    solve_parser.add_argument("--time-limit", type=float)
    solve_parser.add_argument("--node-limit", type=int)
    solve_parser.add_argument("--lower-bound", type=int)
    solve_parser.add_argument("--upper-bound", type=int)
    solve_parser.add_argument(
        "--decision-cap",
        type=int,
        help="fix t to this cap and solve pure feasibility",
    )
    solve_parser.add_argument("--output", type=Path)
    verify_parser = sub.add_parser("verify", help="recompute a saved result")
    verify_parser.add_argument("path", type=Path)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(sys.argv[1:] if argv is None else argv)
    if args.command == "solve":
        return solve(args)
    if args.command == "verify":
        return verify(args.path)
    raise AssertionError(args.command)


if __name__ == "__main__":
    raise SystemExit(main())
