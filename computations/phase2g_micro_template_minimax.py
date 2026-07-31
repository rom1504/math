#!/usr/bin/env python3
"""Exact minimax over the 16 common cross-micro template toggles.

Each input certificate supplies one witness per diagonal at a stated current
micro template.  The script transports every witness energy back to the
untoggled template exactly, pools the witness rounds, and solves

    min_y max_w |E_w(0) - 2 <c_w,y>|,  y in {0,1}^16.

The result controls the finite pooled witness envelope for one common micro
template.  It is not an exact cap result and does not allow a different
template for each diagonal completion.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
from ortools.sat.python import cp_model



def decode_spin(encoded: str) -> np.ndarray:
    bits = int(encoded, 16)
    return np.asarray(
        [1 if bits & (1 << i) else -1 for i in range(56)], dtype=np.int64
    ).reshape(14, 4)


def micro_features(base: np.ndarray, records: list[dict[str, object]]) -> np.ndarray:
    """Vectorized oriented-bundle contributions for all 16 micro positions."""

    spins = np.stack(
        [decode_spin(str(record["spin_bits_little_endian"])) for record in records]
    )
    upper = np.triu(base.astype(np.int64), 1)
    hadamard = np.asarray(
        [[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]],
        dtype=np.int64,
    )
    features = np.empty((len(records), 16), dtype=np.int64)
    for a in range(4):
        for b in range(4):
            features[:, 4 * a + b] = hadamard[a, b] * np.einsum(
                "ri,ij,rj->r",
                spins[:, :, a],
                upper,
                spins[:, :, b],
                optimize=True,
            )
    return features


def load_round(
    path: Path, base: np.ndarray
) -> tuple[np.ndarray, np.ndarray, dict[str, object]]:
    payload = json.loads(path.read_text())
    records = payload["records"]
    if len(records) != 16384 or payload["unresolved_count"] != 0:
        raise AssertionError(f"incomplete round: {path}")
    features = micro_features(base, records)
    current_positions = [
        tuple(int(value) for value in item)
        for item in payload.get("cross_micro_position_toggles", [])
    ]
    current_columns = [4 * a + b for a, b in current_positions]
    saved = np.asarray([int(row["energy"]) for row in records], dtype=np.int64)
    base_energy = saved.copy()
    if current_columns:
        base_energy += 2 * np.sum(features[:, current_columns], axis=1)
    metadata = {
        "path": str(path),
        "record_hash": payload["canonical_record_sha256"],
        "current_positions": [list(item) for item in current_positions],
        "saved_minimum_absolute_energy": int(np.min(np.abs(saved))),
    }
    return base_energy, features.astype(np.int64), metadata


def solve_minimax(
    energies: np.ndarray,
    features: np.ndarray,
    time_limit: float,
    workers: int,
    initial_rows: int,
    rows_per_round: int,
) -> dict[str, object]:
    initial_order = np.argsort(-np.abs(energies), kind="stable")
    active = set(int(row) for row in initial_order[:initial_rows])
    rounds = []
    started = time.time()
    while True:
        remaining = time_limit - (time.time() - started)
        if remaining <= 0:
            return {
                "status": "TIME_LIMIT",
                "rounds": rounds,
                "active_constraint_count": len(active),
            }
        model = cp_model.CpModel()
        toggles = [model.new_bool_var(f"toggle_{column}") for column in range(16)]
        maximum = model.new_int_var(0, 2000, "maximum_absolute_energy")
        for row in sorted(active):
            expression = int(energies[row]) - 2 * sum(
                int(features[row, column]) * toggles[column]
                for column in range(16)
            )
            model.add(expression <= maximum)
            model.add(expression >= -maximum)
        model.minimize(maximum)
        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = remaining
        solver.parameters.num_search_workers = workers
        solver.parameters.random_seed = 20260731
        status = solver.solve(model)
        status_name = solver.status_name(status)
        if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            return {
                "status": status_name,
                "rounds": rounds,
                "active_constraint_count": len(active),
                "objective_lower_bound": solver.best_objective_bound,
            }
        chosen = np.asarray(
            [solver.value(variable) for variable in toggles], dtype=np.int64
        )
        corrected = energies - 2 * (features @ chosen)
        objective = int(solver.value(maximum))
        violations = np.flatnonzero(np.abs(corrected) > objective)
        round_record = {
            "round": len(rounds) + 1,
            "solver_status": status_name,
            "active_constraints": len(active),
            "active_objective": objective,
            "active_objective_lower_bound": int(
                np.ceil(solver.best_objective_bound - 1e-9)
            ),
            "full_objective": int(np.max(np.abs(corrected))),
            "full_violations": len(violations),
            "chosen_columns": np.flatnonzero(chosen).astype(int).tolist(),
        }
        rounds.append(round_record)
        print(round_record, flush=True)
        if len(violations) == 0 and status == cp_model.OPTIMAL:
            return {
                "status": "OPTIMAL",
                "objective": objective,
                "objective_lower_bound": objective,
                "chosen_columns": np.flatnonzero(chosen).astype(int).tolist(),
                "chosen_positions": [
                    [int(column // 4), int(column % 4)]
                    for column in np.flatnonzero(chosen)
                ],
                "minimum_energy": int(np.min(corrected)),
                "maximum_energy": int(np.max(corrected)),
                "active_constraint_count": len(active),
                "rounds": rounds,
            }
        severity = np.abs(corrected[violations]) - objective
        ranked = violations[np.argsort(-severity, kind="stable")]
        active.update(int(row) for row in ranked[:rows_per_round])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("certificates", nargs="+", type=Path)
    parser.add_argument(
        "--source",
        type=Path,
        default=Path("computations/results/heuristic_m14_from_conference.json"),
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--time-limit", type=float, default=600.0)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--initial-rows", type=int, default=256)
    parser.add_argument("--rows-per-round", type=int, default=1024)
    args = parser.parse_args()

    base = np.asarray(json.loads(args.source.read_text())["matrix"], dtype=np.int8)
    loaded = [load_round(path, base) for path in args.certificates]
    energies = np.concatenate([item[0] for item in loaded])
    features = np.concatenate([item[1] for item in loaded], axis=0)
    solved = solve_minimax(
        energies,
        features,
        args.time_limit,
        args.workers,
        args.initial_rows,
        args.rows_per_round,
    )
    payload = {
        "schema": "quadratic-signing-common-micro-template-minimax-v1",
        "classification": (
            "exact finite pooled-witness minimax when status is OPTIMAL; "
            "one common template only, not a cap computation"
        ),
        "rounds": [item[2] for item in loaded],
        "witness_count": len(energies),
        "result": solved,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
