#!/usr/bin/env python3
"""Compute/certify the absolute Boolean cap of one fixed signing with CP-SAT.

For fixed x_0=+1, edge variables encode whether the endpoint spins differ.
The signed quadratic energy is an affine weighted-cut objective.  Two CP-SAT
runs maximize H and -H.  OPTIMAL statuses certify the exact fixed-matrix cap;
FEASIBLE/UNKNOWN runs give rigorous explicit lower witnesses and solver upper
bounds but no exact value.  OR-Tools emits no standalone proof object.
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

from exact_mn_milp import stable_matrix_hash


def load_matrix(path: Path) -> np.ndarray:
    payload = json.loads(path.read_text())
    key = "matrix" if "matrix" in payload else "parent_matrix"
    matrix = np.asarray(payload[key], dtype=np.int8)
    if not np.array_equal(matrix, matrix.T) or np.any(np.diag(matrix)):
        raise ValueError(f"invalid signing matrix in {path}")
    return matrix


def solve_side(
    matrix: np.ndarray,
    objective_sign: int,
    time_limit: float,
    workers: int,
    decision_threshold: int | None,
) -> dict[str, object]:
    n = len(matrix)
    model = cp_model.CpModel()
    spin_bits = [model.new_bool_var(f"y_{i}") for i in range(n)]
    model.add(spin_bits[0] == 0)
    cut_bits = []
    coefficients = []
    constant = 0
    for i, j in combinations(range(n), 2):
        edge = model.new_bool_var(f"cut_{i}_{j}")
        model.add(edge >= spin_bits[i] - spin_bits[j])
        model.add(edge >= spin_bits[j] - spin_bits[i])
        model.add(edge <= spin_bits[i] + spin_bits[j])
        model.add(edge <= 2 - spin_bits[i] - spin_bits[j])
        sign = int(matrix[i, j])
        constant += sign
        cut_bits.append(edge)
        coefficients.append(-2 * sign)
    energy = constant + sum(c * edge for c, edge in zip(coefficients, cut_bits))
    if decision_threshold is None:
        model.maximize(objective_sign * energy)
    else:
        model.add(objective_sign * energy >= decision_threshold)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = workers
    solver.parameters.log_search_progress = True
    solver.parameters.log_to_stdout = True
    started = time.time()
    status = solver.solve(model)
    elapsed = time.time() - started
    status_name = solver.status_name(status)
    result: dict[str, object] = {
        "objective_sign": objective_sign,
        "decision_threshold": decision_threshold,
        "status": status_name,
        "objective": solver.objective_value,
        "best_bound": solver.best_objective_bound,
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
        "wall_time_seconds": solver.wall_time,
        "elapsed_seconds": elapsed,
        "response_stats": solver.response_stats(),
    }
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        spins = np.asarray(
            [1 - 2 * solver.value(bit) for bit in spin_bits], dtype=np.int64
        )
        verified_energy = int(spins @ matrix.astype(np.int64) @ spins // 2)
        verified_objective = objective_sign * verified_energy
        if decision_threshold is None and verified_objective != round(solver.objective_value):
            raise AssertionError((verified_objective, solver.objective_value))
        if decision_threshold is not None and verified_objective < decision_threshold:
            raise AssertionError((verified_objective, decision_threshold))
        result["spins"] = [int(v) for v in spins]
        result["verified_energy"] = verified_energy
    print(
        f"side={objective_sign:+d} status={status_name} "
        f"objective={solver.objective_value} bound={solver.best_objective_bound} "
        f"wall={solver.wall_time:.6f}s",
        flush=True,
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("signing", type=Path)
    parser.add_argument("--time-limit", type=float, default=1800.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--decision-threshold", type=int)
    parser.add_argument("--side", type=int, choices=(-1, 1), action="append")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    matrix = load_matrix(args.signing)
    requested_sides = args.side if args.side else [1, -1]
    sides = [
        solve_side(
            matrix,
            sign,
            args.time_limit,
            args.workers,
            args.decision_threshold,
        )
        for sign in requested_sides
    ]
    certified = args.decision_threshold is None and all(
        side["status"] == "OPTIMAL" for side in sides
    )
    witnessed = [
        abs(int(side["verified_energy"]))
        for side in sides
        if "verified_energy" in side
    ]
    upper_bounds = [float(side["best_bound"]) for side in sides]
    payload = {
        "schema": "quadratic-signing-fixed-cap-cpsat-v1",
        "classification": (
            "solver-certified exact fixed-matrix cap; no standalone proof object"
            if certified
            else (
                (
                    "solver-certified fixed-threshold decision; no standalone proof object"
                    if all(
                        side["status"] in {"OPTIMAL", "FEASIBLE", "INFEASIBLE"}
                        for side in sides
                    )
                    else "inconclusive fixed-threshold solver search"
                )
                if args.decision_threshold is not None
                else "solver bounds and explicit fixed-matrix energy witnesses; not exact"
            )
        ),
        "source": str(args.signing),
        "n": len(matrix),
        "matrix_sha256": stable_matrix_hash(matrix),
        "solver": {
            "ortools_version": ortools.__version__,
            "python_version": platform.python_version(),
            "workers": args.workers,
            "time_limit_seconds_per_side": args.time_limit,
        },
        "decision_threshold": args.decision_threshold,
        "sides": sides,
        "cap_lower_bound": max(witnessed) if witnessed else None,
        "cap_upper_bound": max(upper_bounds) if args.decision_threshold is None else None,
        "certified_cap": max(witnessed) if certified else None,
    }
    print(f"n={len(matrix)} certified={certified} decision={args.decision_threshold}")
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(f"wrote {args.output}")
    valid_statuses = {"OPTIMAL", "FEASIBLE", "INFEASIBLE"}
    return 0 if all(side["status"] in valid_statuses for side in sides) else 2


if __name__ == "__main__":
    raise SystemExit(main())
