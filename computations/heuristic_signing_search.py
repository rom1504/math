#!/usr/bin/env python3
"""Reproducible heuristic search for low-cap quadratic signings.

This is an exploration tool, not an optimality certificate.  It combines
steepest single-edge descent with deterministic-seed random kicks.  Every
reported cap is recomputed exhaustively over all projective Boolean spins.
"""

from __future__ import annotations

import argparse
import json
import time
from itertools import combinations
from pathlib import Path

import numpy as np

from exact_mn_milp import exact_profile, projective_spins, stable_matrix_hash


def search(
    n: int,
    restarts: int,
    kicks: int,
    kick_size: int,
    seed: int,
    target: int | None,
) -> tuple[np.ndarray, dict[str, object]]:
    rng = np.random.default_rng(seed)
    spins = projective_spins(n).astype(np.int16)
    edges = tuple(combinations(range(n), 2))
    products = np.column_stack([spins[:, i] * spins[:, j] for i, j in edges]).astype(
        np.int16
    )
    best_cap = 10**9
    best_signs: np.ndarray | None = None
    trajectory: list[dict[str, int]] = []
    started = time.time()

    for restart in range(restarts):
        signs = rng.choice(np.asarray([-1, 1], dtype=np.int16), size=len(edges))
        energies = products @ signs
        local_kicks = 0
        while local_kicks <= kicks:
            while True:
                candidates = energies[:, None] - 2 * products * signs[None, :]
                caps = np.max(np.abs(candidates), axis=0)
                current = int(np.max(np.abs(energies)))
                edge = int(np.argmin(caps))
                candidate = int(caps[edge])
                if candidate >= current:
                    break
                energies = candidates[:, edge].copy()
                signs[edge] *= -1
            current = int(np.max(np.abs(energies)))
            if current < best_cap:
                best_cap = current
                best_signs = signs.copy()
                trajectory.append(
                    {"restart": restart, "kick": local_kicks, "cap": best_cap}
                )
                print(
                    f"best n={n} cap={best_cap} restart={restart} "
                    f"kick={local_kicks} elapsed={time.time()-started:.3f}s",
                    flush=True,
                )
                if target is not None and best_cap <= target:
                    break
            if local_kicks == kicks:
                break
            chosen = rng.choice(len(edges), size=min(kick_size, len(edges)), replace=False)
            for edge in chosen:
                energies -= 2 * products[:, edge] * signs[edge]
                signs[edge] *= -1
            local_kicks += 1
        if target is not None and best_cap <= target:
            break

    if best_signs is None:
        raise AssertionError("search produced no signing")
    matrix = np.zeros((n, n), dtype=np.int8)
    for (i, j), sign in zip(edges, best_signs):
        matrix[i, j] = matrix[j, i] = int(sign)
    profile = exact_profile(matrix)
    if profile["M"] != best_cap:
        raise AssertionError((profile["M"], best_cap))
    metadata = {
        "schema": "quadratic-signing-heuristic-v1",
        "classification": "heuristic upper bound; exhaustive spin evaluation only",
        "n": n,
        "seed": seed,
        "restarts_requested": restarts,
        "kicks_per_restart": kicks,
        "kick_size": kick_size,
        "target": target,
        "elapsed_seconds": time.time() - started,
        "trajectory": trajectory,
        "matrix": [[int(v) for v in row] for row in matrix],
        "matrix_sha256": stable_matrix_hash(matrix),
        "profile": profile,
    }
    return matrix, metadata


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("n", type=int)
    parser.add_argument("--restarts", type=int, default=100)
    parser.add_argument("--kicks", type=int, default=30)
    parser.add_argument("--kick-size", type=int, default=3)
    parser.add_argument("--seed", type=int, default=20260730)
    parser.add_argument("--target", type=int)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    _, payload = search(
        args.n,
        args.restarts,
        args.kicks,
        args.kick_size,
        args.seed,
        args.target,
    )
    print(
        f"final n={args.n} cap={payload['profile']['M']} "
        f"hash={payload['matrix_sha256']} elapsed={payload['elapsed_seconds']:.3f}s"
    )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
