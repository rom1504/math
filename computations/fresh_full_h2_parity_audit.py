#!/usr/bin/env python3
"""Bounded full-H2 parity audit; lower witnesses and optional CP-SAT bound."""

from __future__ import annotations

import argparse
import json

import numpy as np
from ortools.sat.python import cp_model


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outer-exponent", type=int, default=2)
    parser.add_argument("--restarts", type=int, default=10000)
    parser.add_argument("--seconds", type=float, default=60)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    outer = np.ones((4, 4), dtype=np.int64) - 2 * np.eye(4, dtype=np.int64)
    matrix = np.array([[1, 1], [1, -1]], dtype=np.int64)
    for _ in range(args.outer_exponent):
        matrix = np.kron(outer, matrix)
    n = len(matrix)
    generator = np.random.default_rng(20260905)
    best = -1
    best_u = best_v = None
    for _ in range(args.restarts):
        v = generator.choice((-1, 1), n)
        for iteration in range(100):
            u = np.where(matrix @ v >= 0, 1, -1)
            new_v = np.where(matrix @ u >= 0, 1, -1)
            if np.array_equal(new_v, v):
                break
            v = new_v
        value = int(u @ matrix @ v)
        if value > best:
            best, best_u, best_v = value, u.copy(), v.copy()
    report = {
        "n": n,
        "outer_exponent": args.outer_exponent,
        "bilinear_witness_value": best,
        "bilinear_witness_left": best_u.tolist(),
        "bilinear_witness_right": best_v.tolist(),
        "implied_regularized_lower": best / (2 * (4 ** args.outer_exponent) ** 1.5),
    }
    print(json.dumps(report), flush=True)
    if args.seconds <= 0:
        return
    for orientation in (1, -1):
        model = cp_model.CpModel()
        bits = [model.new_bool_var(f"b_{i}") for i in range(n)]
        model.add(bits[0] == 0)
        objective = []
        offset = int(np.trace(matrix)) // 2
        for i in range(n):
            for j in range(i + 1, n):
                xor = model.new_bool_var(f"xor_{i}_{j}")
                model.add_bool_xor([bits[i], bits[j], xor.Not()])
                coefficient = int(matrix[i, j])
                offset += coefficient
                objective.append(-2 * coefficient * xor)
        energy = sum(objective) + offset
        model.maximize(orientation * energy)
        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = args.seconds
        solver.parameters.num_search_workers = args.workers
        solver.parameters.random_seed = 20260905
        status = solver.solve(model)
        output = {
            "orientation": orientation,
            "status": solver.status_name(status),
            "seconds": solver.wall_time,
            "best_bound": solver.best_objective_bound,
        }
        if status in (cp_model.FEASIBLE, cp_model.OPTIMAL):
            spins = np.array([1 - 2 * solver.value(bit) for bit in bits], dtype=np.int64)
            exact_energy = int(spins @ matrix @ spins) // 2
            assert orientation * exact_energy == round(solver.objective_value)
            output.update({"witness_energy": exact_energy, "witness": spins.tolist()})
        print(json.dumps(output), flush=True)


if __name__ == "__main__":
    main()
