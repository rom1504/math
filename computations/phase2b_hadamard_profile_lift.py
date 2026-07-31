#!/usr/bin/env python3
"""Build low-cap-scale Hadamard lifts of the equal-phi6 collision."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
from pathlib import Path

import numpy as np

from heuristic_fixed_cap_search import search_side
from phase2_restriction_state_audit import class_map, profile


def sylvester(order: int) -> np.ndarray:
    if order < 1 or order & (order - 1):
        raise ValueError("Hadamard order must be a power of two")
    matrix = np.ones((1, 1), dtype=np.int64)
    while len(matrix) < order:
        matrix = np.block([[matrix, matrix], [matrix, -matrix]])
    if not np.array_equal(matrix @ matrix.T, order * np.eye(order, dtype=np.int64)):
        raise AssertionError("Hadamard identity failed")
    return matrix


def lift(base: np.ndarray, hadamard: np.ndarray) -> np.ndarray:
    within = hadamard - np.diag(np.diag(hadamard))
    result = np.kron(base, hadamard) + np.kron(
        np.eye(len(base), dtype=np.int64), within
    )
    if np.any(np.diag(result)) or not np.array_equal(result, result.T):
        raise AssertionError("invalid lift diagonal or symmetry")
    off_diagonal = result[~np.eye(len(result), dtype=bool)]
    if not np.all(np.abs(off_diagonal) == 1):
        raise AssertionError("lift is not a signing")
    return result


def exact_evaluate(executable: Path, matrix: np.ndarray) -> dict[str, object]:
    payload = str(len(matrix)) + "\n" + "\n".join(
        " ".join(map(str, row)) for row in matrix
    ) + "\n"
    result = subprocess.run([str(executable)], input=payload, text=True,
                            capture_output=True, check=True)
    return json.loads(result.stdout)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("collision", type=Path)
    parser.add_argument("--orders", default="2,4,8,16")
    parser.add_argument("--exact-evaluator", type=Path, required=True)
    parser.add_argument("--restarts", type=int, default=400)
    parser.add_argument("--kicks", type=int, default=30)
    parser.add_argument("--kick-size", type=int, default=6)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.collision.read_text())
    bases = [np.asarray(record["matrix"], dtype=np.int64)
             for record in source["records"]]
    labels = {size: class_map(size)[0] for size in (4, 5, 6)}
    records = []
    for order in map(int, args.orders.split(",")):
        hadamard = sylvester(order)
        lifts = [lift(base, hadamard) for base in bases]
        row = {"hadamard_order": order, "lift_order": len(lifts[0]), "lifts": []}
        if order == 2:
            profiles = [
                {str(size): profile(matrix, size, labels[size])
                 for size in (4, 5, 6)} for matrix in lifts
            ]
            if profiles[0] != profiles[1]:
                raise AssertionError("twofold lift profile mismatch")
            row["exact_profile_equality_audit"] = profiles[0]
        for index, matrix in enumerate(lifts):
            rng = np.random.default_rng(args.seed + 1009 * order + index)
            positive = search_side(matrix, args.restarts, args.kicks,
                                   args.kick_size, rng)
            negative = search_side(-matrix, args.restarts, args.kicks,
                                   args.kick_size, rng)
            heuristic_lower = max(int(positive["best_energy"]),
                                  int(negative["best_energy"]))
            eigen_norm = float(np.max(np.abs(np.linalg.eigvalsh(matrix))))
            item = {
                "base_code": source["records"][index]["code"],
                "matrix_sha256": hashlib.sha256(
                    matrix.astype(np.int8).tobytes()).hexdigest(),
                "spectral_norm": eigen_norm,
                "rigorous_spectral_cap_upper": len(matrix) * eigen_norm / 2,
                "heuristic_cap_lower": heuristic_lower,
                "heuristic_normalized_cap_lower": heuristic_lower / len(matrix) ** 1.5,
                "positive_witness": positive,
                "negative_witness": negative,
            }
            if order == 2:
                item["exact_profile"] = exact_evaluate(
                    args.exact_evaluator.resolve(), matrix)
            row["lifts"].append(item)
        records.append(row)
        print(
            f"k={order} heuristic={[r['heuristic_cap_lower'] for r in row['lifts']]} "
            f"spectral={[r['rigorous_spectral_cap_upper'] for r in row['lifts']]}",
            flush=True,
        )
    output = {
        "schema": "quadratic-signing-phase2b-hadamard-profile-lift-v1",
        "classification": (
            "proved explicit signing/profile-preserving construction and rigorous "
            "spectral upper bounds; exact caps only where exact_profile is present; "
            "all other cap lower bounds heuristic with explicit witnesses"
        ),
        "source": str(args.collision),
        "definition": "S_A(H)=A tensor H + I_10 tensor (H-diag(H))",
        "uniform_profile_reason": (
            "oriented switching equivalence of each base support lifts by switching "
            "whole fibers; summing micro-coordinate occupancies is class-dependent only"
        ),
        "uniform_scale_bound": (
            "for symmetric Hadamard H of order k, norm(S)<=norm(A)*sqrt(k)+sqrt(k)+1, "
            "so cap(S)=O((10k)^(3/2))"
        ),
        "heuristic_parameters": {
            "restarts": args.restarts, "kicks": args.kicks,
            "kick_size": args.kick_size, "seed": args.seed,
        },
        "records": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
