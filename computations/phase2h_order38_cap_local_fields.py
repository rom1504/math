#!/usr/bin/env python3
"""Reproducible cap heuristics and local-field audit at order 38.

Compares representatives of all three exactly classified natural cyclic
solution orbits at k=3: one Paley and two graph-inequivalent non-Paley
orbits. Search outcomes are heuristic, but every retained witness energy and
local-field identity is checked by exact integer arithmetic.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np

from audit_two_fiber_boolean_spectral_deficit import (
    gauge_delete_negative_graph,
)
from phase2d_audit_balanced_diagonals import climb_batch
from two_fiber_cyclic_conference import circulant, verify


def matrix_hash(matrix: np.ndarray) -> str:
    return hashlib.sha256(matrix.astype(np.int8).tobytes()).hexdigest()


def matrix_from_sequences(a: list[int], c: list[int]) -> np.ndarray:
    verify(3, a, c)
    matrix_a = circulant(a)
    matrix_c = circulant(c)
    return np.block([[matrix_a, matrix_c], [matrix_c.T, -matrix_a]])


def encode_spin(spin: np.ndarray) -> str:
    bits = sum((int(value) > 0) << i for i, value in enumerate(spin))
    return f"{bits:010x}"


def exact_record(matrix: np.ndarray, spin: np.ndarray, k: int = 3) -> dict[str, object]:
    field = matrix @ spin
    energy = int(spin @ field // 2)
    orientation = 1 if energy >= 0 else -1
    oriented = orientation * spin * field
    penalty = (oriented - (2 * k - 1)) * (oriented - (2 * k + 1))
    deficit = 2 * k * (2 * k * k + 1) - abs(energy)
    if int(np.sum(penalty)) != 8 * k * deficit or np.any(penalty < 0):
        raise AssertionError("local-field penalty identity failed")
    return {
        "spin_positive_bits_little_endian": encode_spin(spin),
        "energy": energy,
        "absolute_energy": abs(energy),
        "deficit_from_114": deficit,
        "oriented_local_field_histogram": {
            str(key): int(value)
            for key, value in sorted(Counter(map(int, oriented)).items())
        },
        "penalty_sum": int(np.sum(penalty)),
    }


def pair_climb(matrix: np.ndarray, spin: np.ndarray) -> tuple[np.ndarray, int, int]:
    """Exact best-improvement ascent over all one- and two-spin flips."""

    spin = spin.copy()
    steps = 0
    while True:
        field = matrix @ spin
        signed = spin * field
        single_delta = -2 * signed
        products = spin[:, None] * spin[None, :]
        pair_delta = (
            -2 * (signed[:, None] + signed[None, :])
            + 4 * matrix * products
        )
        pair_delta[np.tril_indices(len(spin))] = -10**9
        single_index = int(np.argmax(single_delta))
        pair_flat = int(np.argmax(pair_delta))
        left, right = np.unravel_index(pair_flat, pair_delta.shape)
        if int(single_delta[single_index]) <= 0 and int(pair_delta[left, right]) <= 0:
            break
        if int(single_delta[single_index]) >= int(pair_delta[left, right]):
            spin[single_index] *= -1
        else:
            spin[left] *= -1
            spin[right] *= -1
        steps += 1
    energy = int(spin @ matrix @ spin // 2)
    return spin, energy, steps


def search_matrix(
    matrix: np.ndarray,
    seed: int,
    random_starts_per_sign: int,
    batch_size: int,
    pair_starts_per_sign: int,
) -> dict[str, object]:
    rng = np.random.default_rng(seed)
    energy_histogram: Counter[int] = Counter()
    best = -1
    best_records: dict[tuple[int, tuple[tuple[int, int], ...]], dict[str, object]] = {}
    batched_states = 0
    pair_states = 0
    pair_steps = 0
    for sign in (1, -1):
        signed_matrix = sign * matrix
        remaining = random_starts_per_sign
        while remaining:
            count = min(batch_size, remaining)
            initial = rng.choice(
                np.asarray([-1, 1], dtype=np.int64), size=(count, len(matrix))
            )
            spins, energies = climb_batch(signed_matrix, initial)
            for value, multiplicity in zip(*np.unique(energies, return_counts=True)):
                energy_histogram[int(value)] += int(multiplicity)
            batch_best = int(np.max(energies))
            if batch_best >= best:
                for index in np.flatnonzero(energies == batch_best):
                    record = exact_record(matrix, spins[index])
                    magnitude = int(record["absolute_energy"])
                    signature = tuple(
                        (int(key), int(value))
                        for key, value in record[
                            "oriented_local_field_histogram"
                        ].items()
                    )
                    if magnitude > best:
                        best = magnitude
                        best_records.clear()
                    if magnitude == best:
                        best_records[(int(record["energy"]), signature)] = record
            remaining -= count
            batched_states += count

        for _ in range(pair_starts_per_sign):
            initial = rng.choice(
                np.asarray([-1, 1], dtype=np.int64), size=len(matrix)
            )
            spin, energy, steps = pair_climb(signed_matrix, initial)
            pair_steps += steps
            pair_states += 1
            record = exact_record(matrix, spin)
            magnitude = int(record["absolute_energy"])
            signature = tuple(
                (int(key), int(value))
                for key, value in record["oriented_local_field_histogram"].items()
            )
            if magnitude > best:
                best = magnitude
                best_records.clear()
            if magnitude == best:
                best_records[(int(record["energy"]), signature)] = record

    return {
        "seed": seed,
        "random_single_flip_starts_per_sign": random_starts_per_sign,
        "batch_size": batch_size,
        "pair_flip_starts_per_sign": pair_starts_per_sign,
        "single_flip_local_optima_sampled": batched_states,
        "pair_flip_local_optima_sampled": pair_states,
        "pair_climb_total_steps": pair_steps,
        "best_certified_lower_bound": best,
        "single_flip_terminal_energy_histogram": {
            str(key): value for key, value in sorted(energy_histogram.items())
        },
        "distinct_best_energy_field_signatures": list(best_records.values()),
        "classification": (
            "heuristic search; retained witnesses and local-field identities "
            "verified exactly; no cap upper bound"
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--random-starts-per-sign", type=int, default=200000)
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--pair-starts-per-sign", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument(
        "--orbit-source",
        type=Path,
        default=Path("computations/results/two_fiber_cyclic_conference.json"),
    )
    args = parser.parse_args()

    orbit_payload = json.loads(args.orbit_source.read_text())[
        "k3_exhaustive_natural_orbits"
    ]
    if orbit_payload["raw_oriented_solution_count"] != 627:
        raise AssertionError("k=3 exhaustive orbit source changed")
    matrices = {}
    orbit_metadata = {}
    for orbit in orbit_payload["orbits"]:
        counts = orbit["clique_counts"]
        name = (
            f"orbit_K4_{counts['4']}_K5_{counts['5']}"
            + ("_paley" if orbit["is_prime_paley_nonsplit_torus_orbit"] else "_nonpaley")
        )
        matrices[name] = matrix_from_sequences(orbit["a"], orbit["c"])
        orbit_metadata[name] = {
            "natural_orbit_raw_size": orbit["natural_orbit_raw_size"],
            "clique_counts": counts,
            "is_prime_paley_nonsplit_torus_orbit": orbit[
                "is_prime_paley_nonsplit_torus_orbit"
            ],
            "a": orbit["a"],
            "c": orbit["c"],
        }
    results = {}
    for offset, (name, matrix) in enumerate(matrices.items()):
        if not np.array_equal(matrix @ matrix, 37 * np.eye(38, dtype=np.int64)):
            raise AssertionError(f"conference identity failed for {name}")
        results[name] = {
            "orbit": orbit_metadata[name],
            "matrix_sha256": matrix_hash(matrix),
            "root0_negative_graph_degree_histogram": {
                str(key): int(value)
                for key, value in sorted(
                    Counter(
                        map(int, np.sum(gauge_delete_negative_graph(matrix), axis=1))
                    ).items()
                )
            },
            "search": search_matrix(
                matrix,
                args.seed + offset,
                args.random_starts_per_sign,
                args.batch_size,
                args.pair_starts_per_sign,
            ),
        }
        print(
            f"{name}: best={results[name]['search']['best_certified_lower_bound']}",
            flush=True,
        )

    payload = {
        "schema": "quadratic-signing-order38-local-field-comparison-v1",
        "classification": (
            "reproducible heuristic cap comparison with exact witness checks; "
            "no exact cap or scaling claim"
        ),
        "universal_exact_upper_bound": 113,
        "results": results,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
