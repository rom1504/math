#!/usr/bin/env python3
"""Probe low-variance vertices of the exact fractional bridge body.

This is an exploratory, reproducible LP experiment.  It enumerates every
projective child spin pair, builds the polytope (A.10) from
artifacts/second_phase_independent_abstraction.md, and optimizes independent
random linear objectives.  Each optimum is a vertex candidate; the script
records the smallest rounding variance sum(1-c_ij^2) found.  It does *not*
certify the globally smallest variance, since that is a convex-quadratic
maximization problem.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


def spins(n: int) -> np.ndarray:
    vals = np.arange(1 << (n - 1), dtype=np.uint64)[:, None]
    bits = ((vals >> np.arange(n - 1, dtype=np.uint64)) & 1).astype(np.int8)
    return np.concatenate((np.ones((len(vals), 1), dtype=np.int8), 1 - 2 * bits), axis=1)


def energy(matrix: np.ndarray, states: np.ndarray) -> np.ndarray:
    states64 = states.astype(np.int64, copy=False)
    matrix64 = matrix.astype(np.int64, copy=False)
    return np.einsum("bi,ij,bj->b", states64, matrix64, states64, optimize=True) // 2


def load_matrix(path: Path) -> np.ndarray:
    data = json.loads(path.read_text())
    return np.asarray(data["matrix"], dtype=np.int8)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--child-a", type=Path, required=True)
    parser.add_argument("--child-b", type=Path, required=True)
    parser.add_argument("--sign-b", type=int, choices=(-1, 1), default=1)
    parser.add_argument("--target", type=float, required=True)
    parser.add_argument("--trials", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    a = load_matrix(args.child_a)
    b = args.sign_b * load_matrix(args.child_b)
    m, n = len(a), len(b)
    xs, ys = spins(m), spins(n)
    ha, hb = energy(a, xs), energy(b, ys)

    # One row vector vec(x y^T) for every projective pair.
    vectors = np.einsum("ai,bj->abij", xs, ys, optimize=True).reshape(-1, m * n).astype(float)
    margins = (args.target - np.abs(ha[:, None] + hb[None, :])).reshape(-1).astype(float)
    if np.min(margins) < -1e-9:
        raise SystemExit("fractional body empty already at C_0=0: a margin is negative")
    aub = np.concatenate((vectors, -vectors), axis=0)
    bub = np.concatenate((margins, margins), axis=0)

    rng = np.random.default_rng(args.seed)
    best: dict[str, object] | None = None
    statuses: dict[str, int] = {}
    for trial in range(args.trials):
        objective = rng.normal(size=m * n)
        result = linprog(-objective, A_ub=aub, b_ub=bub, bounds=(-1.0, 1.0), method="highs")
        statuses[str(result.status)] = statuses.get(str(result.status), 0) + 1
        if not result.success:
            continue
        c = result.x
        variance = float(np.sum(1.0 - c * c))
        fractional = np.flatnonzero(np.abs(c) < 1.0 - 1e-7)
        residual = aub @ c - bub
        active = np.flatnonzero(np.abs(residual) <= 1e-7)
        if len(fractional):
            active_rank = int(np.linalg.matrix_rank(aub[np.ix_(active, fractional)], tol=1e-8))
        else:
            active_rank = 0
        record = {
            "trial": trial,
            "variance": variance,
            "variance_fraction": variance / (m * n),
            "fractional_entries": int(len(fractional)),
            "active_state_inequalities": int(len(active)),
            "active_rank_on_fractional_entries": active_rank,
            "max_constraint_residual": float(np.max(residual)),
            "bridge": c.reshape(m, n).tolist(),
        }
        if best is None or variance < float(best["variance"]):
            best = record

    payload = {
        "schema": "quadratic-signing-fractional-bridge-variance-probe-v1",
        "classification": "heuristic LP vertex search; all spin-pair constraints are exact",
        "child_a": str(args.child_a),
        "child_b": str(args.child_b),
        "sign_b": args.sign_b,
        "orders": [m, n],
        "target": args.target,
        "projective_state_pairs": int(len(vectors)),
        "bridge_variables": m * n,
        "linear_inequalities": int(len(aub)),
        "trials": args.trials,
        "seed": args.seed,
        "solver_status_counts": statuses,
        "best": best,
        "warning": "random LP objectives do not certify the global minimum of variance",
    }
    encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(encoded)
    print(encoded, end="")
    print("sha256", hashlib.sha256(encoded.encode()).hexdigest())


if __name__ == "__main__":
    main()
