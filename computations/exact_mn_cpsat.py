#!/usr/bin/env python3
"""CP-SAT decision/optimization model for the exact quadratic-signing cap.

This is an independent backend for the rooted, symmetry-reduced model from
`exact_mn_milp.py`, with an additional valid degree ordering inside the two
remaining neighbor groups.  Its most useful mode is a fixed-cap infeasibility
test:

    exact_mn_cpsat.py 11 --decision-cap 15

An INFEASIBLE solver status certifies computationally that no signing reaches
the proposed cap within the encoded symmetry-complete search space.  OR-Tools
does not emit a standalone proof certificate here, so the classification is
solver-certified computation rather than a formal proof.
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


def build_model(n: int, decision_cap: int | None, root_degree: int | None = None) -> tuple[
    cp_model.CpModel,
    tuple[tuple[int, int], ...],
    list[cp_model.IntVar],
    cp_model.IntVar | None,
]:
    if n < 3:
        raise ValueError("n must be at least 3")
    model = cp_model.CpModel()
    edges = tuple(combinations(range(1, n), 2))
    edge_index = {edge: k for k, edge in enumerate(edges)}
    z = [model.new_bool_var(f"z_{i}_{j}") for i, j in edges]

    # The same complement symmetry permits choosing the rooted internal graph
    # with no more than half of its possible negative edges.
    model.add(sum(z) <= len(z) // 2)

    # Symmetry-complete rooted constraints, identical to the MILP backend.
    model.add(z[edge_index[(1, 2)]] == 0)
    for j in range(2, n - 1):
        model.add(z[edge_index[(1, j)]] <= z[edge_index[(1, j + 1)]])
    incident_1 = [z[edge_index[(1, j)]] for j in range(2, n)]
    negative_degrees: dict[int, cp_model.LinearExpr] = {1: sum(incident_1)}
    for i in range(2, n):
        incident_i = []
        for j in range(1, n):
            if i == j:
                continue
            edge = (j, i) if j < i else (i, j)
            incident_i.append(z[edge_index[edge]])
        negative_degrees[i] = sum(incident_i)
        model.add(negative_degrees[1] <= negative_degrees[i])

    # Complementing every rooted internal edge is induced by globally
    # negating the signing and re-gauging the root.  For a graph on n-1
    # vertices, either it or its complement has minimum degree at most
    # floor((n-2)/2).  Choose that representative.  Optional equality permits
    # an exact disjoint case split across the remaining root degrees.
    model.add(negative_degrees[1] <= (n - 2) // 2)
    if root_degree is not None:
        if not 0 <= root_degree <= (n - 2) // 2:
            raise ValueError("root degree is outside the symmetry-complete range")
        model.add(negative_degrees[1] == root_degree)

    # The vertices after 1 have already been split into its positive and
    # negative neighbors.  Permutations inside either group remain free, so
    # sort their negative degrees.  The two enforcement alternatives cover
    # equal adjacent incident signs; the 0-to-1 boundary needs no ordering.
    for j in range(2, n - 1):
        left = z[edge_index[(1, j)]]
        right = z[edge_index[(1, j + 1)]]
        model.add(negative_degrees[j] <= negative_degrees[j + 1]).only_enforce_if(
            [left.Not(), right.Not()]
        )
        model.add(negative_degrees[j] <= negative_degrees[j + 1]).only_enforce_if(
            [left, right]
        )

    if decision_cap is None:
        cap = model.new_int_var(0, n * (n - 1) // 2, "cap")
        model.minimize(cap)
    else:
        cap = None

    for mask in range(1 << (n - 1)):
        tail = [1 - 2 * ((mask >> j) & 1) for j in range(n - 1)]
        constant = sum(tail)
        coefficients = []
        for i, j in edges:
            product = tail[i - 1] * tail[j - 1]
            constant += product
            coefficients.append(-2 * product)
        energy = constant + sum(c * variable for c, variable in zip(coefficients, z))
        bound = decision_cap if decision_cap is not None else cap
        model.add(energy <= bound)
        model.add(energy >= -bound)
    return model, edges, z, cap


def matrix_from_values(
    n: int,
    edges: tuple[tuple[int, int], ...],
    values: list[int],
) -> np.ndarray:
    matrix = np.zeros((n, n), dtype=np.int8)
    matrix[0, 1:] = matrix[1:, 0] = 1
    for (i, j), value in zip(edges, values):
        sign = 1 - 2 * int(value)
        matrix[i, j] = matrix[j, i] = sign
    return matrix


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("n", type=int)
    parser.add_argument("--decision-cap", type=int)
    parser.add_argument("--root-degree", type=int)
    parser.add_argument("--time-limit", type=float, default=1800.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    model, edges, variables, cap = build_model(
        args.n, args.decision_cap, args.root_degree
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
        f"status={status_name} objective={solver.objective_value} "
        f"best_bound={solver.best_objective_bound} "
        f"conflicts={solver.num_conflicts} branches={solver.num_branches} "
        f"wall={solver.wall_time:.6f}s elapsed={elapsed:.6f}s",
        flush=True,
    )

    payload: dict[str, object] = {
        "schema": "quadratic-signing-exact-cpsat-v1",
        "classification": "solver-certified computation; no standalone proof object",
        "n": args.n,
        "normalization": "M_n=max_x |sum_{i<j} a_ij x_i x_j|",
        "mode": "decision" if args.decision_cap is not None else "optimization",
        "decision_cap": args.decision_cap,
        "model": {
            "root_gauge": True,
            "basic_permutation_and_complement_symmetry": True,
            "root_negative_degree_case": args.root_degree,
            "internal_binary_variables": len(edges),
            "projective_spin_constraints": 2 * (1 << (args.n - 1)),
        },
        "solver": {
            "ortools_version": ortools.__version__,
            "python_version": platform.python_version(),
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
        values = [solver.value(variable) for variable in variables]
        matrix = matrix_from_values(args.n, edges, values)
        profile = exact_profile(matrix)
        if args.decision_cap is not None and profile["M"] > args.decision_cap:
            raise AssertionError((profile["M"], args.decision_cap))
        if cap is not None and profile["M"] > round(solver.value(cap)):
            raise AssertionError((profile["M"], solver.value(cap)))
        payload["matrix"] = [[int(v) for v in row] for row in matrix]
        payload["matrix_sha256"] = stable_matrix_hash(matrix)
        payload["profile"] = profile
        print(
            f"verified profile M={profile['M']} P={profile['P']} Q={profile['Q']} "
            f"sha256={payload['matrix_sha256']}",
            flush=True,
        )

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(f"wrote {args.output}", flush=True)
    return 0 if status in (cp_model.OPTIMAL, cp_model.FEASIBLE, cp_model.INFEASIBLE) else 2


if __name__ == "__main__":
    raise SystemExit(main())
