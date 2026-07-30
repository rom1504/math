#!/usr/bin/env python3
"""Optimize a bridge block between two fixed quadratic signings.

For child matrices A and B, this solves exactly

  min_C max_{x,y} (|H_A(x)+H_B(y)| + |x^T C y|),

where C is an m-by-n sign matrix.  By independently flipping the global sign
of one child spin, this is the absolute quadratic cap of the assembled block
signing.  The model therefore tests the state-dependent composition problem,
not merely the standalone bipartite norm of C.

The result is solver-certified computation.  The returned assembled signing
is exhaustively evaluated and hashed, but OR-Tools does not emit a standalone
formal optimality certificate.
"""

from __future__ import annotations

import argparse
import json
import time
from itertools import product
from pathlib import Path

import numpy as np
import ortools
from ortools.sat.python import cp_model

from exact_mn_milp import exact_profile, projective_spins, stable_matrix_hash


def load_matrix(path: Path, class_index: int | None = None) -> np.ndarray:
    payload = json.loads(path.read_text())
    matrix_data = (
        payload["matrix"]
        if class_index is None
        else payload["classes"][class_index]["representative_matrix"]
    )
    matrix = np.asarray(matrix_data, dtype=np.int8)
    if not np.array_equal(matrix, matrix.T) or np.any(np.diag(matrix)):
        raise ValueError(f"invalid signing matrix in {path}")
    return matrix


def one_copy_energies(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    spins = projective_spins(len(matrix)).astype(np.int64)
    energies = np.einsum("bi,ij,bj->b", spins, matrix.astype(np.int64), spins) // 2
    return spins, energies


def build_model(
    a: np.ndarray,
    b: np.ndarray,
    sign_b: int,
    decision_cap: int | None,
) -> tuple[cp_model.CpModel, list[cp_model.IntVar], cp_model.IntVar]:
    m, n = len(a), len(b)
    spins_a, energies_a = one_copy_energies(a)
    spins_b, energies_b = one_copy_energies(sign_b * b)
    model = cp_model.CpModel()
    variables = [model.new_bool_var(f"z_{i}_{j}") for i in range(m) for j in range(n)]
    # C and -C have the same objective.
    model.add(variables[0] == 0)
    total_edges = (m + n) * (m + n - 1) // 2
    if decision_cap is None:
        cap = model.new_int_var(0, total_edges, "cap")
        model.minimize(cap)
    else:
        cap = model.new_constant(decision_cap)

    for ia, ib in product(range(len(spins_a)), range(len(spins_b))):
        x = spins_a[ia]
        y = spins_b[ib]
        constant = int(x.sum() * y.sum())
        coefficients = []
        for i in range(m):
            for j in range(n):
                coefficients.append(-2 * int(x[i] * y[j]))
        cross = constant + sum(c * z for c, z in zip(coefficients, variables))
        internal = abs(int(energies_a[ia] + energies_b[ib]))
        model.add(cross + internal <= cap)
        model.add(-cross + internal <= cap)
    return model, variables, cap


def assemble(a: np.ndarray, b: np.ndarray, sign_b: int, bridge: np.ndarray) -> np.ndarray:
    m, n = len(a), len(b)
    matrix = np.zeros((m + n, m + n), dtype=np.int8)
    matrix[:m, :m] = a
    matrix[m:, m:] = sign_b * b
    matrix[:m, m:] = bridge
    matrix[m:, :m] = bridge.T
    return matrix


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("child_a", type=Path)
    parser.add_argument("child_b", type=Path)
    parser.add_argument("--child-a-class", type=int)
    parser.add_argument("--child-b-class", type=int)
    parser.add_argument("--sign-b", type=int, choices=(-1, 1), default=1)
    parser.add_argument("--decision-cap", type=int)
    parser.add_argument("--time-limit", type=float, default=600.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    a = load_matrix(args.child_a, args.child_a_class)
    b = load_matrix(args.child_b, args.child_b_class)
    model, variables, cap = build_model(a, b, args.sign_b, args.decision_cap)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = args.time_limit
    solver.parameters.num_search_workers = args.workers
    solver.parameters.log_search_progress = True
    solver.parameters.log_to_stdout = True
    started = time.time()
    status = solver.solve(model)
    elapsed = time.time() - started
    status_name = solver.status_name(status)
    print(
        f"status={status_name} objective={solver.objective_value} "
        f"bound={solver.best_objective_bound} branches={solver.num_branches} "
        f"wall={solver.wall_time:.6f}s",
        flush=True,
    )

    payload: dict[str, object] = {
        "schema": "quadratic-signing-bridge-cpsat-v1",
        "classification": "solver-certified computation; no standalone proof object",
        "child_a": str(args.child_a),
        "child_b": str(args.child_b),
        "child_a_class": args.child_a_class,
        "child_b_class": args.child_b_class,
        "sign_b": args.sign_b,
        "orders": [len(a), len(b)],
        "mode": "decision" if args.decision_cap is not None else "optimization",
        "decision_cap": args.decision_cap,
        "solver": {
            "ortools_version": ortools.__version__,
            "status": status_name,
            "objective": solver.objective_value,
            "best_bound": solver.best_objective_bound,
            "conflicts": solver.num_conflicts,
            "branches": solver.num_branches,
            "wall_time_seconds": solver.wall_time,
            "elapsed_seconds": elapsed,
            "workers": args.workers,
            "time_limit_seconds": args.time_limit,
            "response_stats": solver.response_stats(),
        },
    }

    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        m, n = len(a), len(b)
        values = np.asarray([solver.value(z) for z in variables], dtype=np.int8)
        bridge = (1 - 2 * values).reshape(m, n)
        parent = assemble(a, b, args.sign_b, bridge)
        profile_a = exact_profile(a)
        profile_b = exact_profile(b)
        profile_parent = exact_profile(parent)
        if args.decision_cap is None and profile_parent["M"] != round(solver.value(cap)):
            raise AssertionError((profile_parent["M"], solver.value(cap)))
        if args.decision_cap is not None and profile_parent["M"] > args.decision_cap:
            raise AssertionError((profile_parent["M"], args.decision_cap))
        lhs = float(profile_parent["M"] ** (2.0 / 3.0))
        rhs = float(profile_a["M"] ** (2.0 / 3.0) + profile_b["M"] ** (2.0 / 3.0))
        payload.update(
            {
                "bridge": [[int(v) for v in row] for row in bridge],
                "parent_matrix": [[int(v) for v in row] for row in parent],
                "parent_matrix_sha256": stable_matrix_hash(parent),
                "child_profiles": [profile_a, profile_b],
                "parent_profile": profile_parent,
                "two_thirds_power": {
                    "parent": lhs,
                    "child_sum": rhs,
                    "defect": lhs - rhs,
                },
            }
        )
        print(
            f"verified bridge {m}+{n}: parent_M={profile_parent['M']} "
            f"two_thirds_defect={lhs-rhs:+.12f} "
            f"hash={payload['parent_matrix_sha256']}",
            flush=True,
        )

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(f"wrote {args.output}", flush=True)
    return 0 if status in (cp_model.OPTIMAL, cp_model.FEASIBLE, cp_model.INFEASIBLE) else 2


if __name__ == "__main__":
    raise SystemExit(main())
