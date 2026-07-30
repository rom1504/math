#!/usr/bin/env python3
"""Reproducible coordinate-ascent lower bounds for one fixed signing cap.

This is an exploratory maximization tool, not an optimality certificate.  It
performs exact single-spin ascent for both energy signs, with reproducible
random restarts and kicks.  Every saved witness is checked by direct integer
matrix multiplication, so it is a rigorous lower bound on the cap of the fixed
matrix even though the search for the maximum is heuristic.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np

from exact_mn_milp import stable_matrix_hash


def climb(matrix: np.ndarray, spins: np.ndarray) -> tuple[np.ndarray, int, int]:
    spins = spins.copy()
    field = matrix @ spins
    energy = int(spins @ field // 2)
    flips = 0
    while True:
        signed_fields = spins * field
        vertex = int(np.argmin(signed_fields))
        if int(signed_fields[vertex]) >= 0:
            break
        old_spin = int(spins[vertex])
        energy -= 2 * int(signed_fields[vertex])
        field -= 2 * old_spin * matrix[:, vertex]
        spins[vertex] = -old_spin
        flips += 1
    verified = int(spins @ matrix @ spins // 2)
    if verified != energy:
        raise AssertionError((verified, energy))
    return spins, energy, flips


def search_side(
    matrix: np.ndarray,
    restarts: int,
    kicks: int,
    kick_size: int,
    rng: np.random.Generator,
) -> dict[str, object]:
    n = len(matrix)
    best_energy = -10**30
    best_spins = None
    total_climb_flips = 0
    improvements = []
    for restart in range(restarts):
        spins = rng.choice(np.asarray([-1, 1], dtype=np.int64), size=n)
        for kick in range(kicks + 1):
            spins, energy, flips = climb(matrix, spins)
            total_climb_flips += flips
            if energy > best_energy:
                best_energy = energy
                best_spins = spins.copy()
                improvements.append(
                    {"restart": restart, "kick": kick, "energy": best_energy}
                )
                print(
                    f"improvement restart={restart} kick={kick} energy={energy}",
                    flush=True,
                )
            if kick < kicks:
                chosen = rng.choice(n, size=min(kick_size, n), replace=False)
                spins = spins.copy()
                spins[chosen] *= -1
    if best_spins is None:
        raise AssertionError("empty search")
    return {
        "best_energy": int(best_energy),
        "spins": [int(value) for value in best_spins],
        "total_climb_flips": total_climb_flips,
        "improvements": improvements,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument(
        "--matrix-key",
        choices=("matrix", "parent_matrix", "conference_matrix"),
    )
    parser.add_argument("--restarts", type=int, default=200)
    parser.add_argument("--kicks", type=int, default=20)
    parser.add_argument("--kick-size", type=int, default=4)
    parser.add_argument("--seed", type=int, default=20260730)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = json.loads(args.source.read_text())
    key = args.matrix_key
    if key is None:
        key = "matrix" if "matrix" in payload else "parent_matrix"
    matrix = np.asarray(payload[key], dtype=np.int64)
    if not np.array_equal(matrix, matrix.T) or np.any(np.diag(matrix)):
        raise ValueError("selected matrix is not a symmetric zero-diagonal signing")
    rng = np.random.default_rng(args.seed)
    started = time.time()
    positive = search_side(
        matrix, args.restarts, args.kicks, args.kick_size, rng
    )
    negative = search_side(
        -matrix, args.restarts, args.kicks, args.kick_size, rng
    )
    positive_energy = int(positive["best_energy"])
    negative_energy = int(negative["best_energy"])
    cap_lower = max(positive_energy, negative_energy)
    output = {
        "schema": "quadratic-signing-fixed-cap-heuristic-v1",
        "classification": (
            "heuristic search with explicit exactly verified energy witnesses; "
            "rigorous fixed-matrix cap lower bound only"
        ),
        "source": str(args.source),
        "matrix_key": key,
        "n": len(matrix),
        "matrix_sha256": stable_matrix_hash(matrix),
        "seed": args.seed,
        "restarts": args.restarts,
        "kicks": args.kicks,
        "kick_size": args.kick_size,
        "positive": positive,
        "negative": negative,
        "cap_lower_bound": cap_lower,
        "elapsed_seconds": time.time() - started,
    }
    print(
        f"n={len(matrix)} positive={positive_energy} negative={negative_energy} "
        f"cap_lower={cap_lower} elapsed={output['elapsed_seconds']:.3f}s"
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
        print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
