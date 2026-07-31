#!/usr/bin/env python3
"""Heuristically maximize squared norm of a feasible fractional bridge.

For fixed child signings A,B and target T, bridge entries c_ij in [-1,1]
obey the exact linear constraints

  |H_A(x)+H_B(y)+x^T C y| <= T

for every projective parent spin.  The rounding variance is
V=sum(1-c_ij^2), so minimizing V is convex maximization over this polytope.
We search exposed vertices with random objectives and monotone tangent ascent.
The feasible points and their measured V are reproducible numerical evidence;
global optimality of V is not claimed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import time
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


def all_spins(n: int, fix_first: bool) -> np.ndarray:
    free = n - 1 if fix_first else n
    codes = np.arange(1 << free, dtype=np.uint64)
    bits = ((codes[:, None] >> np.arange(free, dtype=np.uint64)) & 1).astype(np.int8)
    spins = 1 - 2 * bits
    if fix_first:
        spins = np.column_stack([np.ones(len(spins), dtype=np.int8), spins])
    return spins


def cap(matrix: np.ndarray) -> int:
    spins = all_spins(len(matrix), True).astype(np.int64)
    energies = np.einsum("bi,ij,bj->b", spins, matrix, spins) // 2
    return int(np.max(np.abs(energies)))


def nested_value(payload: object, key: str) -> object:
    value = payload
    for component in key.split("."):
        value = value[int(component)] if isinstance(value, list) else value[component]
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("left", type=Path)
    parser.add_argument("right", type=Path)
    parser.add_argument("--left-key", default="matrix")
    parser.add_argument("--right-key", default="matrix")
    parser.add_argument("--restarts", type=int, default=12)
    parser.add_argument("--iterations", type=int, default=8)
    parser.add_argument("--target-slack", type=float, default=0.0)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    left_payload = json.loads(args.left.read_text())
    right_payload = json.loads(args.right.read_text())
    a = np.asarray(nested_value(left_payload, args.left_key), dtype=np.int64)
    b = np.asarray(nested_value(right_payload, args.right_key), dtype=np.int64)
    m, n = len(a), len(b)
    ma, mb = cap(a), cap(b)
    ideal_target = (ma ** (2 / 3) + mb ** (2 / 3)) ** 1.5
    target = ideal_target + args.target_slack

    xs = all_spins(m, True).astype(np.float64)
    ys = all_spins(n, False).astype(np.float64)
    internal_x = np.einsum("bi,ij,bj->b", xs, a, xs) / 2
    internal_y = np.einsum("bi,ij,bj->b", ys, b, ys) / 2
    internal = (internal_x[:, None] + internal_y[None, :]).reshape(-1)
    products = np.einsum("ai,bj->abij", xs, ys).reshape(-1, m * n)
    aub = np.vstack([products, -products])
    bub = np.concatenate([target - internal, target + internal])
    rng = np.random.default_rng(args.seed)
    best = None
    trajectories = []
    started = time.monotonic()
    for restart in range(args.restarts):
        objective = rng.normal(size=m * n)
        row = []
        solution = None
        for iteration in range(args.iterations + 1):
            result = linprog(-objective, A_ub=aub, b_ub=bub,
                             bounds=[(-1, 1)] * (m * n), method="highs")
            if not result.success:
                raise RuntimeError((result.status, result.message))
            solution = result.x
            squared = float(solution @ solution)
            variance = m * n - squared
            row.append({"iteration": iteration, "squared_norm": squared,
                        "variance": variance})
            objective = 2 * solution
        assert solution is not None
        if best is None or float(solution @ solution) > best[0]:
            best = (float(solution @ solution), solution.copy(), restart)
        trajectories.append(row)
        print(f"restart={restart} V={m*n-float(solution@solution):.9f}", flush=True)
    assert best is not None
    squared, solution, best_restart = best
    parent_energies = internal + products @ solution
    numerical_cap = float(np.max(np.abs(parent_energies)))
    rounded = np.where(solution >= 0, 1, -1).reshape(m, n)
    rounded_parent = np.block([[a, rounded], [rounded.T, b]])
    rounded_cap = cap(rounded_parent)
    output = {
        "schema": "quadratic-signing-phase2-fractional-bridge-variance-v1",
        "classification": (
            "reproducible heuristic convex-maximization over an exactly "
            "enumerated fractional bridge polytope; no global V optimum claimed"
        ),
        "left": str(args.left), "right": str(args.right),
        "orders": [m, n], "child_caps": [ma, mb],
        "ideal_power_target": ideal_target,
        "target_slack": args.target_slack,
        "fractional_target": target,
        "projective_parent_constraints": len(internal),
        "fractional_bridge": solution.reshape(m, n).tolist(),
        "fractional_bridge_sha256": hashlib.sha256(
            solution.astype(np.float64).tobytes()).hexdigest(),
        "fractional_numerical_cap": numerical_cap,
        "maximum_constraint_violation": max(0.0, numerical_cap - target),
        "squared_norm": squared,
        "rounding_variance": m * n - squared,
        "rounding_variance_over_bridge_edges": (m * n - squared) / (m * n),
        "coordinates_with_abs_at_least_1_minus_1e-8":
            int(np.count_nonzero(np.abs(solution) >= 1 - 1e-8)),
        "rounded_bridge_cap_exact": rounded_cap,
        "best_restart": best_restart,
        "seed": args.seed, "restarts": args.restarts,
        "tangent_iterations": args.iterations,
        "trajectories": trajectories,
        "elapsed_seconds": time.monotonic() - started,
        "scipy_method": "HiGHS via scipy.optimize.linprog",
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(f"best V={m*n-squared:.9f} target={target:.9f} rounded_cap={rounded_cap}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
