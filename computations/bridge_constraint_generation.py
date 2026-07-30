#!/usr/bin/env python3
"""Exact fixed-child bridge search by CP-SAT constraint generation."""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import ortools
from ortools.sat.python import cp_model

from conference_prime_square import evaluate, stable_hash


def spins_from_code(code: int, n: int) -> np.ndarray:
    spins = np.ones(n, dtype=np.int64)
    for vertex in range(1, n):
        if code & (1 << (vertex - 1)):
            spins[vertex] = -1
    return spins


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--matrix-key", default="conference_matrix")
    parser.add_argument("--left-vertices", required=True)
    parser.add_argument("--target-cap", type=int, required=True)
    parser.add_argument("--evaluator", type=Path, required=True)
    parser.add_argument("--max-iterations", type=int, default=500)
    parser.add_argument("--solve-time", type=float, default=10.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument(
        "--bundle-extremizers",
        action="store_true",
        help="separate every global maximizer and minimizer of each candidate",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--resume",
        type=Path,
        help="preload separated_gray_codes from an earlier compatible result",
    )
    args = parser.parse_args()
    payload = json.loads(args.input.read_text())
    source = np.asarray(payload[args.matrix_key], dtype=np.int8)
    left_vertices = tuple(int(item) for item in args.left_vertices.split(","))
    right_vertices = tuple(i for i in range(len(source)) if i not in left_vertices)
    left = source[np.ix_(left_vertices, left_vertices)].astype(np.int64)
    right = source[np.ix_(right_vertices, right_vertices)].astype(np.int64)
    original_bridge = source[np.ix_(left_vertices, right_vertices)].astype(np.int64)
    m, n = len(left), len(right)
    total_order = m + n

    model = cp_model.CpModel()
    variables = [[model.new_bool_var(f"z_{i}_{j}") for j in range(n)] for i in range(m)]
    flip_terms = []
    for i in range(m):
        for j in range(n):
            flip_terms.append(variables[i][j] if original_bridge[i, j] == 1 else 1 - variables[i][j])
    model.minimize(sum(flip_terms))
    seen_codes = set()
    iterations = []
    started = time.monotonic()

    def add_spin_constraint(code: int) -> None:
        if code in seen_codes:
            return
        seen_codes.add(code)
        spin = spins_from_code(code, total_order)
        x, y = spin[:m], spin[m:]
        internal = int((x @ left @ x + y @ right @ y) // 2)
        products = np.outer(x, y)
        constant = internal + int(products.sum())
        expression = constant + sum(
            -2 * int(products[i, j]) * variables[i][j]
            for i in range(m)
            for j in range(n)
        )
        model.add(expression <= args.target_cap)
        model.add(expression >= -args.target_cap)

    # Begin with both extremal states of the inherited conference bridge.
    initial_parent = np.block([[left, original_bridge], [original_bridge.T, right]])
    initial_profile = evaluate(
        args.evaluator.resolve(), initial_parent, args.bundle_extremizers
    )
    initial_codes = [
        int(initial_profile["argmax_gray"]),
        int(initial_profile["argmin_gray"]),
    ]
    if args.bundle_extremizers:
        initial_codes = [
            *map(int, initial_profile["maximizer_gray_codes"]),
            *map(int, initial_profile["minimizer_gray_codes"]),
        ]
    for code in initial_codes:
        add_spin_constraint(code)
    if args.resume is not None:
        resume_payload = json.loads(args.resume.read_text())
        if resume_payload["left_vertices"] != list(left_vertices):
            raise ValueError("resume file has different fixed children")
        if resume_payload["target_cap"] != args.target_cap:
            raise ValueError("resume file has a different target cap")
        for code in resume_payload["separated_gray_codes"]:
            add_spin_constraint(int(code))

    final_parent = None
    final_profile = None
    final_status = "ITERATION_LIMIT"
    for iteration in range(1, args.max_iterations + 1):
        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = args.solve_time
        solver.parameters.num_search_workers = args.workers
        status = solver.solve(model)
        status_name = solver.status_name(status)
        row = {
            "iteration": iteration,
            "constraints": len(seen_codes),
            "solver_status": status_name,
            "solver_wall_time_seconds": solver.wall_time,
        }
        if status == cp_model.INFEASIBLE:
            final_status = "INFEASIBLE"
            iterations.append(row)
            print(f"iteration={iteration} constraints={len(seen_codes)} INFEASIBLE", flush=True)
            break
        if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            final_status = status_name
            iterations.append(row)
            print(f"iteration={iteration} constraints={len(seen_codes)} {status_name}", flush=True)
            break
        bridge = np.asarray(
            [[1 - 2 * solver.value(variables[i][j]) for j in range(n)] for i in range(m)],
            dtype=np.int64,
        )
        parent = np.block([[left, bridge], [bridge.T, right]])
        profile = evaluate(args.evaluator.resolve(), parent, args.bundle_extremizers)
        row.update(
            {
                "candidate_cap": profile["cap"],
                "candidate_min_energy": profile["min_energy"],
                "candidate_max_energy": profile["max_energy"],
                "bridge_flips_from_seed": int(np.count_nonzero(bridge != original_bridge)),
                "matrix_sha256": profile["matrix_sha256"],
            }
        )
        iterations.append(row)
        print(
            f"iteration={iteration} constraints={len(seen_codes)} cap={profile['cap']} "
            f"range=[{profile['min_energy']},{profile['max_energy']}] "
            f"flips={row['bridge_flips_from_seed']}",
            flush=True,
        )
        if int(profile["cap"]) <= args.target_cap:
            final_status = "FEASIBLE"
            final_parent = parent
            final_profile = profile
            break
        before = len(seen_codes)
        codes = [int(profile["argmax_gray"]), int(profile["argmin_gray"])]
        if args.bundle_extremizers:
            codes = [
                *map(int, profile["maximizer_gray_codes"]),
                *map(int, profile["minimizer_gray_codes"]),
            ]
        for code in codes:
            add_spin_constraint(code)
        row["new_extremal_state_constraints"] = len(seen_codes) - before
        if len(seen_codes) == before:
            raise AssertionError("violating extrema were already constrained")

    output = {
        "schema": "quadratic-signing-bridge-constraint-generation-v1",
        "classification": (
            "solver-certified fixed-child bridge computation with exact exhaustive candidate separation; no standalone proof object"
        ),
        "input": str(args.input),
        "orders": [m, n],
        "left_vertices": list(left_vertices),
        "right_vertices": list(right_vertices),
        "target_cap": args.target_cap,
        "initial_profile": initial_profile,
        "initial_matrix_sha256": stable_hash(initial_parent),
        "ortools_version": ortools.__version__,
        "max_iterations": args.max_iterations,
        "solve_time_seconds_per_iteration": args.solve_time,
        "workers": args.workers,
        "bundle_extremizers": args.bundle_extremizers,
        "elapsed_seconds": time.monotonic() - started,
        "final_status": final_status,
        "iterations": iterations,
        "separated_projective_states": len(seen_codes),
        "separated_gray_codes": sorted(seen_codes),
        "evaluator_source": "computations/exact_fixed_signing_gray.cpp",
    }
    if final_parent is not None and final_profile is not None:
        output["parent_matrix"] = [[int(value) for value in row] for row in final_parent]
        output["parent_matrix_sha256"] = stable_hash(final_parent)
        output["parent_profile"] = final_profile
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(f"final_status={final_status} wrote {args.output}")
    return 0 if final_status in {"FEASIBLE", "INFEASIBLE"} else 2


if __name__ == "__main__":
    raise SystemExit(main())
