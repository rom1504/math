#!/usr/bin/env python3
"""Certify caps of Sylvester-Hadamard replicas of exact small minimizers.

For every signed-permutation/global-sign class recorded by
``enumerate_minimizer_orbits.py``, this script forms ``H_k \otimes A`` and
solves both signs of its Boolean quadratic maximum by CP-SAT.  An OPTIMAL
status on both sides certifies the exact fixed-core cap (the solver emits no
standalone proof object).  The intended use is a finite falsification test for
the selectable-replica inequality

    min_{A in Argmin M_r, k in K} [M(H_k tensor A)-k^(3/2) M_r] <= 0.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import time
from itertools import combinations
from pathlib import Path

import numpy as np
import ortools
from ortools.sat.python import cp_model

from exact_mn_milp import stable_matrix_hash


ROOT = Path(__file__).resolve().parents[1]


def sylvester(order: int) -> np.ndarray:
    if order < 1 or order & (order - 1):
        raise ValueError("Sylvester order must be a positive power of two")
    matrix = np.asarray([[1]], dtype=np.int8)
    seed = np.asarray([[1, 1], [1, -1]], dtype=np.int8)
    while len(matrix) < order:
        matrix = np.kron(matrix, seed).astype(np.int8)
    return matrix


def solve_side(
    matrix: np.ndarray,
    objective_sign: int,
    time_limit: float,
    workers: int,
) -> dict[str, object]:
    n = len(matrix)
    model = cp_model.CpModel()
    spin_bits = [model.new_bool_var(f"y_{i}") for i in range(n)]
    model.add(spin_bits[0] == 0)
    constant = 0
    terms = []
    for i, j in combinations(range(n), 2):
        edge = model.new_bool_var(f"cut_{i}_{j}")
        model.add(edge >= spin_bits[i] - spin_bits[j])
        model.add(edge >= spin_bits[j] - spin_bits[i])
        model.add(edge <= spin_bits[i] + spin_bits[j])
        model.add(edge <= 2 - spin_bits[i] - spin_bits[j])
        coefficient = int(matrix[i, j])
        constant += coefficient
        terms.append(-2 * coefficient * edge)
    energy = constant + sum(terms)
    model.maximize(objective_sign * energy)

    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = time_limit
    solver.parameters.num_search_workers = workers
    solver.parameters.log_search_progress = False
    started = time.time()
    status = solver.solve(model)
    elapsed = time.time() - started
    answer: dict[str, object] = {
        "objective_sign": objective_sign,
        "status": solver.status_name(status),
        "objective": solver.objective_value,
        "best_bound": solver.best_objective_bound,
        "conflicts": solver.num_conflicts,
        "branches": solver.num_branches,
        "wall_time_seconds": solver.wall_time,
        "elapsed_seconds": elapsed,
    }
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        spins = np.asarray(
            [1 - 2 * solver.value(bit) for bit in spin_bits], dtype=np.int64
        )
        verified_energy = int(spins @ matrix.astype(np.int64) @ spins // 2)
        if objective_sign * verified_energy != round(solver.objective_value):
            raise AssertionError((verified_energy, solver.objective_value))
        answer["verified_energy"] = verified_energy
        answer["spins"] = [int(value) for value in spins]
    return answer


def load_classes(order: int) -> tuple[int, list[dict[str, object]]]:
    source = ROOT / "computations" / "results" / f"m{order}_minimizer_orbits.json"
    payload = json.loads(source.read_text())
    return int(payload["target_cap"]), payload["classes"]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orders", type=int, nargs="+", default=list(range(3, 9)))
    parser.add_argument("--ks", type=int, nargs="+", default=[2, 4])
    parser.add_argument("--time-limit", type=float, default=1800.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "computations" / "results" / "hadamard_replica_cap_selection.json",
    )
    args = parser.parse_args()

    records: list[dict[str, object]] = []
    for order in args.orders:
        child_cap, classes = load_classes(order)
        for orbit in classes:
            child = np.asarray(orbit["representative_matrix"], dtype=np.int8)
            for replica_order in args.ks:
                core = np.kron(sylvester(replica_order), child).astype(np.int8)
                sides = [
                    solve_side(core, sign, args.time_limit, args.workers)
                    for sign in (1, -1)
                ]
                certified = all(side["status"] == "OPTIMAL" for side in sides)
                witnessed = max(
                    abs(int(side["verified_energy"]))
                    for side in sides
                    if "verified_energy" in side
                )
                cap = witnessed if certified else None
                threshold = replica_order ** 1.5 * child_cap
                record = {
                    "child_order": order,
                    "child_class": int(orbit["class"]),
                    "child_cap": child_cap,
                    "child_matrix_sha256": stable_matrix_hash(child),
                    "replica_order": replica_order,
                    "core_order": len(core),
                    "core_matrix_sha256": stable_matrix_hash(core),
                    "sides": sides,
                    "certified_core_cap": cap,
                    "core_cap_lower_bound": witnessed,
                    "ideal_core_cap": threshold,
                    "core_excess": None if cap is None else cap - threshold,
                    "b_defect_core_only": (
                        None
                        if cap is None
                        else cap ** (2 / 3) - replica_order * child_cap ** (2 / 3)
                    ),
                }
                records.append(record)
                print(
                    f"r={order} class={orbit['class']} k={replica_order} "
                    f"N={len(core)} cap={cap or ('>=' + str(witnessed))} "
                    f"excess={record['core_excess']}",
                    flush=True,
                )

    summary: list[dict[str, object]] = []
    for order in args.orders:
        candidates = [record for record in records if record["child_order"] == order]
        certified_candidates = [
            record for record in candidates if record["certified_core_cap"] is not None
        ]
        best = min(certified_candidates, key=lambda record: float(record["core_excess"]))
        summary.append(
            {
                "child_order": order,
                "best_child_class": best["child_class"],
                "best_replica_order": best["replica_order"],
                "best_core_cap": best["certified_core_cap"],
                "best_core_excess": best["core_excess"],
                "selectable_nonpositive": float(best["core_excess"]) <= 0,
            }
        )

    output = {
        "schema": "exact-small-minimizer-hadamard-replica-cap-selection-v1",
        "classification": "CP-SAT-certified exact fixed-core caps; no standalone proof object",
        "normalization": "M(A)=max_x |sum_(i<j) a_ij x_i x_j|",
        "solver": {
            "ortools_version": ortools.__version__,
            "python_version": platform.python_version(),
            "workers": args.workers,
            "time_limit_seconds_per_side": args.time_limit,
        },
        "records": records,
        "selection_summary": summary,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
