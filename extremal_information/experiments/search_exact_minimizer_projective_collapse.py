#!/usr/bin/env python3
"""Search for an exact-cap signing with a projectively localized top shell.

This is a finite falsification probe, not an asymptotic theorem and not a
standalone proof-certificate generator.  We gauge one positive ground spin to
the all-one word.  Every projective spin farther than ``vertex_radius`` from
that word is forced below the absolute ``shell_deficit`` shell.

Example:

    .venv/bin/python extremal_information/experiments/\
      search_exact_minimizer_projective_collapse.py \
      10 11 --vertex-radius 1 --shell-deficit 0 --time-limit 60
"""

from __future__ import annotations

import argparse
import itertools
import json
import platform
import time
from pathlib import Path

import numpy as np
import ortools
from ortools.sat.python import cp_model


def projective_spins(n: int):
    return [(1,) + tail for tail in itertools.product((-1, 1), repeat=n - 1)]


def build_model(n: int, cap: int, vertex_radius: int, shell_deficit: int):
    edges = tuple(itertools.combinations(range(n), 2))
    model = cp_model.CpModel()
    negative = [model.new_bool_var(f"z_{i}_{j}") for i, j in edges]

    # The declared ground is all-one and positive.  This fixes switching and
    # global coefficient-negation gauges without losing any candidate.
    model.add(sum(negative) == (len(edges) - cap) // 2)

    degrees = []
    for i in range(n):
        degrees.append(sum(z for z, (u, v) in zip(negative, edges) if i in (u, v)))
    # Permuting vertices preserves all constraints, so degree ordering is a
    # valid (modest) symmetry reduction.
    for i in range(n - 1):
        model.add(degrees[i] <= degrees[i + 1])

    far_bound = cap - shell_deficit - 2
    if far_bound < 0:
        raise ValueError("shell deficit leaves no parity-compatible far energy")

    for spin in projective_spins(n):
        products = [spin[i] * spin[j] for i, j in edges]
        energy = sum(products) - 2 * sum(p * z for p, z in zip(products, negative))
        model.add(energy <= cap)
        model.add(energy >= -cap)
        weight = sum(value < 0 for value in spin)
        distance = min(weight, n - weight)
        if distance > vertex_radius:
            model.add(energy <= far_bound)
            model.add(energy >= -far_bound)

    return model, edges, negative


def matrix_from_solution(n: int, edges, negative, solver):
    matrix = np.zeros((n, n), dtype=np.int8)
    for (i, j), variable in zip(edges, negative):
        value = 1 - 2 * solver.value(variable)
        matrix[i, j] = matrix[j, i] = value
    return matrix


def exact_profile(matrix: np.ndarray, shell_deficit: int):
    n = len(matrix)
    rows = []
    for spin in projective_spins(n):
        vector = np.asarray(spin, dtype=np.int64)
        energy = int(vector @ matrix.astype(np.int64) @ vector // 2)
        rows.append((spin, energy))
    cap = max(abs(energy) for _, energy in rows)
    shell_floor = cap - shell_deficit
    shell_rows = [
        (spin, energy) for spin, energy in rows if abs(energy) >= shell_floor
    ]
    return {
        "cap": cap,
        "shell_deficit": shell_deficit,
        "shell_absolute_energy_floor": shell_floor,
        "positive_ground_count": sum(energy == cap for _, energy in rows),
        "negative_ground_count": sum(energy == -cap for _, energy in rows),
        "absolute_ground_vertex_radii": sorted(
            min(sum(v < 0 for v in spin), n - sum(v < 0 for v in spin))
            for spin, energy in rows if abs(energy) == cap
        ),
        "absolute_shell_count": len(shell_rows),
        "absolute_shell_vertex_radii": sorted(
            min(sum(v < 0 for v in spin), n - sum(v < 0 for v in spin))
            for spin, _ in shell_rows
        ),
        "absolute_energy_histogram": {
            str(value): sum(abs(energy) == value for _, energy in rows)
            for value in sorted({abs(energy) for _, energy in rows})
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("n", type=int)
    parser.add_argument("cap", type=int)
    parser.add_argument("--vertex-radius", type=int, required=True)
    parser.add_argument("--shell-deficit", type=int, default=0)
    parser.add_argument("--time-limit", type=float, default=300.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument(
        "--log-search", action="store_true", help="emit verbose CP-SAT progress"
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    edge_count = args.n * (args.n - 1) // 2
    if (edge_count - args.cap) % 2:
        raise ValueError("cap parity must equal edge-count parity")
    if args.shell_deficit < 0 or args.shell_deficit % 2:
        raise ValueError("shell deficit must be a nonnegative even integer")

    model, edges, negative = build_model(
        args.n, args.cap, args.vertex_radius, args.shell_deficit
    )
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = args.time_limit
    solver.parameters.num_search_workers = args.workers
    solver.parameters.log_search_progress = args.log_search
    solver.parameters.log_to_stdout = args.log_search
    started = time.time()
    status = solver.solve(model)
    elapsed = time.time() - started

    payload = {
        "schema": "exact-minimizer-projective-collapse-v1",
        "classification": "solver-certified finite computation; no standalone proof object",
        "n": args.n,
        "cap": args.cap,
        "vertex_radius": args.vertex_radius,
        "shell_deficit": args.shell_deficit,
        "far_absolute_energy_bound": args.cap - args.shell_deficit - 2,
        "solver": {
            "ortools_version": ortools.__version__,
            "python_version": platform.python_version(),
            "status": solver.status_name(status),
            "wall_time_seconds": solver.wall_time,
            "elapsed_seconds": elapsed,
            "conflicts": solver.num_conflicts,
            "branches": solver.num_branches,
            "workers": args.workers,
            "time_limit_seconds": args.time_limit,
        },
    }
    if status in (cp_model.FEASIBLE, cp_model.OPTIMAL):
        matrix = matrix_from_solution(args.n, edges, negative, solver)
        profile = exact_profile(matrix, args.shell_deficit)
        if profile["cap"] != args.cap:
            raise AssertionError(profile)
        if max(profile["absolute_ground_vertex_radii"]) > args.vertex_radius:
            raise AssertionError(profile)
        if max(profile["absolute_shell_vertex_radii"]) > args.vertex_radius:
            raise AssertionError(profile)
        payload["matrix"] = matrix.astype(int).tolist()
        payload["verified_profile"] = profile

    print(json.dumps(payload, indent=2, sort_keys=True))
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return 0 if status in (cp_model.FEASIBLE, cp_model.OPTIMAL, cp_model.INFEASIBLE) else 2


if __name__ == "__main__":
    raise SystemExit(main())
