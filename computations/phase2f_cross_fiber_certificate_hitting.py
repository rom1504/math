#!/usr/bin/env python3
"""Find minimum common cross-fiber corrections killing saved bad witnesses.

The input contains one exactly verified order-56 witness for every diagonal
completion of a fixed order-14 seed.  For a common set of cross-fiber edge
toggles, a witness energy changes by

    E_D' = E_D - 2 sum_e c[D,e] y_e,

where c[D,e] is the signed contribution of edge e to the saved witness.
This script minimizes either the number of individual toggled edges or the
number of whole 4-by-4 macro-edge blocks subject to |E_D'| <= 208 for every
saved witness.

It uses exact CP-SAT constraint generation.  Once an optimum for the active
rows satisfies every omitted row, that optimum is also globally optimal:
the active model supplies the lower bound and the full check supplies
feasibility.  The result concerns the finite saved witness set, not the true
caps of the corrected signings.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from ortools.sat.python import cp_model


def sylvester4() -> np.ndarray:
    return np.asarray(
        [
            [1, 1, 1, 1],
            [1, -1, 1, -1],
            [1, 1, -1, -1],
            [1, -1, -1, 1],
        ],
        dtype=np.int8,
    )


def decode_spin(encoded: str) -> np.ndarray:
    bits = int(encoded, 16)
    return np.asarray(
        [1 if bits & (1 << i) else -1 for i in range(56)], dtype=np.int8
    ).reshape(14, 4)


def edge_features(
    base: np.ndarray, records: list[dict[str, object]]
) -> tuple[np.ndarray, list[tuple[int, int, int, int]]]:
    edges = [
        (i, a, j, b)
        for i in range(14)
        for j in range(i + 1, 14)
        for a in range(4)
        for b in range(4)
    ]
    features = np.empty((len(records), len(edges)), dtype=np.int8)
    hadamard = sylvester4()
    for row, record in enumerate(records):
        spin = decode_spin(str(record["spin_bits_little_endian"]))
        features[row] = np.asarray(
            [
                base[i, j] * hadamard[a, b] * spin[i, a] * spin[j, b]
                for i, a, j, b in edges
            ],
            dtype=np.int8,
        )
    return features, edges


def block_features(
    edge_matrix: np.ndarray,
    edges: list[tuple[int, int, int, int]],
) -> tuple[np.ndarray, list[tuple[int, int]]]:
    blocks = [(i, j) for i in range(14) for j in range(i + 1, 14)]
    columns = []
    for i, j in blocks:
        selected = [
            column
            for column, (ii, _a, jj, _b) in enumerate(edges)
            if ii == i and jj == j
        ]
        if len(selected) != 16:
            raise AssertionError((i, j, len(selected)))
        columns.append(np.sum(edge_matrix[:, selected], axis=1))
    return np.stack(columns, axis=1).astype(np.int16), blocks


def micro_position_features(
    edge_matrix: np.ndarray,
    edges: list[tuple[int, int, int, int]],
) -> tuple[np.ndarray, list[tuple[int, int]]]:
    """Bundle one oriented micro position across all 91 macro edges."""

    positions = [(a, b) for a in range(4) for b in range(4)]
    columns = []
    for a, b in positions:
        selected = [
            column
            for column, (_i, aa, _j, bb) in enumerate(edges)
            if aa == a and bb == b
        ]
        if len(selected) != 91:
            raise AssertionError((a, b, len(selected)))
        columns.append(np.sum(edge_matrix[:, selected], axis=1))
    return np.stack(columns, axis=1).astype(np.int16), positions


def solve_by_constraint_generation(
    features: np.ndarray,
    energies: np.ndarray,
    target: int,
    time_limit: float,
    workers: int,
    initial_rows: int,
    rows_per_round: int,
) -> dict[str, object]:
    lower_rhs = (energies - target) // 2
    upper_rhs = (energies + target) // 2
    if not np.array_equal(2 * lower_rhs, energies - target):
        raise AssertionError("energy parity mismatch")
    if not np.array_equal(2 * upper_rhs, energies + target):
        raise AssertionError("energy parity mismatch")
    order = np.argsort(-energies, kind="stable")
    active = set(int(row) for row in order[:initial_rows])
    rounds = []
    started = time.time()
    incumbent = None
    objective_bound = None
    while True:
        elapsed = time.time() - started
        remaining = time_limit - elapsed
        if remaining <= 0:
            return {
                "status": "TIME_LIMIT",
                "rounds": rounds,
                "active_constraint_count": len(active),
                "incumbent": incumbent,
                "objective_lower_bound": objective_bound,
            }
        model = cp_model.CpModel()
        variables = [
            model.new_bool_var(f"toggle_{column}")
            for column in range(features.shape[1])
        ]
        for row in sorted(active):
            expression = sum(
                int(features[row, column]) * variables[column]
                for column in range(features.shape[1])
            )
            model.add(expression >= int(lower_rhs[row]))
            model.add(expression <= int(upper_rhs[row]))
        model.minimize(sum(variables))
        solver = cp_model.CpSolver()
        solver.parameters.max_time_in_seconds = remaining
        solver.parameters.num_search_workers = workers
        solver.parameters.random_seed = 20260731
        status = solver.solve(model)
        status_name = solver.status_name(status)
        objective_bound = int(np.ceil(solver.best_objective_bound - 1e-9))
        if status not in (cp_model.OPTIMAL, cp_model.FEASIBLE):
            return {
                "status": status_name,
                "rounds": rounds,
                "active_constraint_count": len(active),
                "objective_lower_bound": objective_bound,
            }
        chosen = np.asarray(
            [solver.value(variable) for variable in variables], dtype=np.int8
        )
        incumbent = np.flatnonzero(chosen).astype(int).tolist()
        sums = features @ chosen.astype(features.dtype)
        corrected = energies - 2 * sums.astype(np.int64)
        violated = np.flatnonzero(np.abs(corrected) > target)
        round_record = {
            "round": len(rounds) + 1,
            "solver_status": status_name,
            "active_constraints": len(active),
            "objective": len(incumbent),
            "objective_lower_bound": objective_bound,
            "full_violations": len(violated),
            "maximum_corrected_energy": int(np.max(corrected)),
            "minimum_corrected_energy": int(np.min(corrected)),
        }
        rounds.append(round_record)
        print(round_record, flush=True)
        if len(violated) == 0 and status == cp_model.OPTIMAL:
            return {
                "status": "OPTIMAL",
                "rounds": rounds,
                "active_constraint_count": len(active),
                "objective": len(incumbent),
                "objective_lower_bound": objective_bound,
                "chosen_columns": incumbent,
                "minimum_corrected_energy": int(np.min(corrected)),
                "maximum_corrected_energy": int(np.max(corrected)),
                "corrected_energy_sha256": hashlib.sha256(
                    corrected.astype("<i8").tobytes()
                ).hexdigest(),
            }
        # A merely feasible active solution supplies no exact lower bound.
        # Continue only after adding the most severe omitted violations.
        severity = np.abs(corrected[violated]) - target
        ranked = violated[np.argsort(-severity, kind="stable")]
        before = len(active)
        active.update(int(row) for row in ranked[:rows_per_round])
        if len(active) == before:
            return {
                "status": "INTERNAL_STALL",
                "rounds": rounds,
                "incumbent": incumbent,
            }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--certificate",
        type=Path,
        default=Path("computations/results/phase2e_all_diagonal_family_audit.json"),
    )
    parser.add_argument(
        "--source",
        type=Path,
        default=Path("computations/results/heuristic_m14_from_conference.json"),
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--representation",
        choices=("edges", "blocks", "micro", "both"),
        default="both",
    )
    parser.add_argument("--target", type=int, default=208)
    parser.add_argument("--time-limit", type=float, default=1800.0)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--initial-rows", type=int, default=64)
    parser.add_argument("--rows-per-round", type=int, default=128)
    args = parser.parse_args()

    certificate = json.loads(args.certificate.read_text())
    source = json.loads(args.source.read_text())
    records = certificate["records"]
    if len(records) != 1 << 14 or certificate["unresolved_count"] != 0:
        raise AssertionError("expected complete all-diagonal certificate")
    base = np.asarray(source["matrix"], dtype=np.int8)
    energies = np.asarray([int(row["energy"]) for row in records], dtype=np.int64)
    edge_matrix, edges = edge_features(base, records)
    results = {}
    representations = (
        ("edges", "blocks", "micro")
        if args.representation == "both"
        else (args.representation,)
    )
    for representation in representations:
        if representation == "edges":
            features = edge_matrix
            objects = edges
            cost_per_object = 1
        elif representation == "blocks":
            features, objects = block_features(edge_matrix, edges)
            cost_per_object = 16
        else:
            features, objects = micro_position_features(edge_matrix, edges)
            cost_per_object = 91
        print(
            f"solving representation={representation} variables={len(objects)}",
            flush=True,
        )
        solved = solve_by_constraint_generation(
            features,
            energies,
            args.target,
            args.time_limit,
            args.workers,
            args.initial_rows,
            args.rows_per_round,
        )
        if "chosen_columns" in solved:
            solved["chosen_objects"] = [
                list(objects[column]) for column in solved["chosen_columns"]
            ]
            solved["modified_edge_count"] = (
                int(solved["objective"]) * cost_per_object
            )
            solved["modified_cross_edge_fraction"] = (
                solved["modified_edge_count"] / len(edges)
            )
        results[representation] = solved

    payload = {
        "schema": "quadratic-signing-cross-fiber-certificate-hitting-v1",
        "classification": (
            "exact CP-SAT optimum when status is OPTIMAL, checked against all "
            "saved witnesses; this kills only the certified finite witness set "
            "and is not a cap upper bound"
        ),
        "certificate": str(args.certificate),
        "certificate_canonical_record_sha256": certificate[
            "canonical_record_sha256"
        ],
        "source": str(args.source),
        "target_absolute_energy": args.target,
        "witness_count": len(records),
        "cross_fiber_edge_count": len(edges),
        "results": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
