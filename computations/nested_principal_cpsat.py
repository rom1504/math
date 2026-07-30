#!/usr/bin/env python3
"""CP-SAT feasibility for a low-cap parent with a low-cap principal child.

The first `child_order` vertices induce the distinguished child.  Vertex 0 is
put in the usual switching gauge, and global negation followed by regauging
fixes edge (1,2) positive.  Permutations are otherwise left available because
the distinguished child makes the full-model ordering constraints invalid.

An OPTIMAL (satisfaction) result supplies an exhaustively verified witness.
INFEASIBLE is a solver-certified finite non-nesting result; OR-Tools emits no
standalone proof object.
"""

from __future__ import annotations

import argparse
import json
import platform
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import ortools
from ortools.sat.python import cp_model

from exact_mn_milp import exact_profile, stable_matrix_hash


def add_cap_constraints(
    model: cp_model.CpModel,
    order: int,
    cap: int,
    edge_index: dict[tuple[int, int], int],
    variables: list[cp_model.IntVar],
) -> None:
    internal_edges = tuple(combinations(range(1, order), 2))
    for mask in range(1 << (order - 1)):
        tail = [1 - 2 * ((mask >> j) & 1) for j in range(order - 1)]
        constant = sum(tail)
        terms = []
        for i, j in internal_edges:
            product = tail[i - 1] * tail[j - 1]
            constant += product
            terms.append(-2 * product * variables[edge_index[(i, j)]])
        energy = constant + sum(terms)
        model.add(energy <= cap)
        model.add(energy >= -cap)


def matrix_from_values(
    order: int,
    edges: tuple[tuple[int, int], ...],
    values: list[int],
) -> np.ndarray:
    matrix = np.zeros((order, order), dtype=np.int8)
    matrix[0, 1:] = matrix[1:, 0] = 1
    for (i, j), value in zip(edges, values):
        matrix[i, j] = matrix[j, i] = 1 - 2 * int(value)
    return matrix


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("parent_order", type=int)
    parser.add_argument("child_order", type=int)
    parser.add_argument("--parent-cap", type=int, required=True)
    parser.add_argument("--child-cap", type=int, required=True)
    parser.add_argument("--time-limit", type=float, default=1800.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if not 3 <= args.child_order < args.parent_order:
        raise ValueError("require 3 <= child_order < parent_order")

    edges = tuple(combinations(range(1, args.parent_order), 2))
    edge_index = {edge: index for index, edge in enumerate(edges)}
    model = cp_model.CpModel()
    variables = [model.new_bool_var(f"z_{i}_{j}") for i, j in edges]
    model.add(variables[edge_index[(1, 2)]] == 0)
    add_cap_constraints(
        model, args.parent_order, args.parent_cap, edge_index, variables
    )
    add_cap_constraints(
        model, args.child_order, args.child_cap, edge_index, variables
    )

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
        f"status={status_name} conflicts={solver.num_conflicts} "
        f"branches={solver.num_branches} wall={solver.wall_time:.6f}s",
        flush=True,
    )
    payload: dict[str, object] = {
        "schema": "quadratic-signing-nested-principal-cpsat-v1",
        "classification": "solver-certified finite computation; no standalone proof object",
        "parent_order": args.parent_order,
        "child_order": args.child_order,
        "parent_cap": args.parent_cap,
        "child_cap": args.child_cap,
        "model": {
            "root_gauge": True,
            "global_complement_symmetry": True,
            "binary_variables": len(variables),
            "parent_projective_constraints": 2 * (1 << (args.parent_order - 1)),
            "child_projective_constraints": 2 * (1 << (args.child_order - 1)),
        },
        "solver": {
            "ortools_version": ortools.__version__,
            "python_version": platform.python_version(),
            "status": status_name,
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
        matrix = matrix_from_values(
            args.parent_order,
            edges,
            [solver.value(variable) for variable in variables],
        )
        parent_profile = exact_profile(matrix)
        child_profile = exact_profile(matrix[: args.child_order, : args.child_order])
        if parent_profile["M"] > args.parent_cap or child_profile["M"] > args.child_cap:
            raise AssertionError((parent_profile["M"], child_profile["M"]))
        payload.update(
            {
                "matrix": [[int(v) for v in row] for row in matrix],
                "matrix_sha256": stable_matrix_hash(matrix),
                "parent_profile": parent_profile,
                "child_profile": child_profile,
            }
        )
        print(
            f"verified parent_M={parent_profile['M']} child_M={child_profile['M']} "
            f"hash={payload['matrix_sha256']}",
            flush=True,
        )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(f"wrote {args.output}", flush=True)
    return 0 if status in (cp_model.OPTIMAL, cp_model.FEASIBLE, cp_model.INFEASIBLE) else 2


if __name__ == "__main__":
    raise SystemExit(main())
