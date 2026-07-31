#!/usr/bin/env python3
"""Audit diagonal completions in the order-14, order-56 Hadamard lift.

For the saved order-14 conference signing ``A`` and every diagonal
``D in {+1,-1}^14`` (balanced by default), form

    S_D = A tensor H_4 + D tensor (H_4 - diag(H_4)).

The search is heuristic, but every successful record contains a Boolean
witness whose energy is recomputed by exact integer arithmetic.  Thus a
complete run with no holdouts is an exhaustive finite *family* certificate:
it does not certify the cap of an individual matrix, only the claimed common
lower bound.  Unresolved diagonals are always reported as open holdouts.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import time
from pathlib import Path

import numpy as np


def sylvester4() -> np.ndarray:
    return np.asarray(
        [
            [1, 1, 1, 1],
            [1, -1, 1, -1],
            [1, 1, -1, -1],
            [1, -1, -1, 1],
        ],
        dtype=np.int64,
    )


def matrix_hash(matrix: np.ndarray) -> str:
    return hashlib.sha256(matrix.astype(np.int8).tobytes()).hexdigest()


def candidate_diagonals(family: str) -> list[tuple[int, np.ndarray]]:
    records = []
    if family == "balanced":
        for positive in itertools.combinations(range(14), 7):
            diagonal = -np.ones(14, dtype=np.int64)
            diagonal[list(positive)] = 1
            mask = sum(1 << vertex for vertex in positive)
            records.append((mask, diagonal))
    elif family == "all":
        for mask in range(1 << 14):
            diagonal = np.asarray(
                [1 if mask & (1 << i) else -1 for i in range(14)],
                dtype=np.int64,
            )
            records.append((mask, diagonal))
    else:
        raise ValueError(f"unknown diagonal family: {family}")
    return records


def lift(
    base: np.ndarray,
    diagonal: np.ndarray,
    cross_micro_position_toggles: tuple[tuple[int, int], ...] = (),
) -> np.ndarray:
    hadamard = sylvester4()
    off_hadamard = hadamard - np.diag(np.diag(hadamard))
    matrix = np.kron(base, hadamard) + np.kron(
        np.diag(diagonal), off_hadamard
    )
    for a, b in cross_micro_position_toggles:
        if not (0 <= a < 4 and 0 <= b < 4):
            raise ValueError((a, b))
        for i in range(14):
            for j in range(i + 1, 14):
                u = 4 * i + a
                v = 4 * j + b
                matrix[u, v] *= -1
                matrix[v, u] *= -1
    if matrix.shape != (56, 56):
        raise AssertionError(matrix.shape)
    if not np.array_equal(matrix, matrix.T) or np.any(np.diag(matrix)):
        raise AssertionError("lift is not symmetric and zero diagonal")
    if not np.all(np.abs(matrix[~np.eye(56, dtype=bool)]) == 1):
        raise AssertionError("lift is not a signing")
    return matrix


def exact_energy(matrix: np.ndarray, spin: np.ndarray) -> int:
    quadratic = int(spin @ matrix @ spin)
    if quadratic % 2:
        raise AssertionError(quadratic)
    return quadratic // 2


def encode_spin(spin: np.ndarray) -> str:
    bits = sum((int(value) > 0) << i for i, value in enumerate(spin))
    return f"{bits:014x}"


def decode_spin(encoded: str) -> np.ndarray:
    bits = int(encoded, 16)
    return np.asarray(
        [1 if bits & (1 << i) else -1 for i in range(56)], dtype=np.int64
    )


def climb_batch(
    matrix: np.ndarray, spins: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Run exact steepest single-spin ascent on every row of ``spins``."""

    spins = spins.copy()
    fields = spins @ matrix
    energies = np.sum(spins * fields, axis=1) // 2
    rows = np.arange(len(spins))
    while True:
        signed_fields = spins * fields
        vertices = np.argmin(signed_fields, axis=1)
        chosen = signed_fields[rows, vertices]
        active = np.flatnonzero(chosen < 0)
        if len(active) == 0:
            break
        active_vertices = vertices[active]
        old = spins[active, active_vertices].copy()
        energies[active] -= 2 * chosen[active]
        fields[active] -= (
            2 * old[:, None] * matrix[active_vertices, :]
        )
        spins[active, active_vertices] = -old
    verified = np.einsum("bi,ij,bj->b", spins, matrix, spins) // 2
    if not np.array_equal(verified, energies):
        raise AssertionError("batched coordinate-ascent bookkeeping failed")
    return spins, energies


def search_witness(
    matrix: np.ndarray,
    rng: np.random.Generator,
    threshold: int,
    batches: int,
    batch_size: int,
) -> tuple[np.ndarray | None, int | None, int]:
    """Search both energy signs and return an exactly checked witness."""

    best_spin = None
    best_energy = -1
    trials = 0
    for _ in range(batches):
        initial = rng.choice(
            np.asarray([-1, 1], dtype=np.int64), size=(batch_size, 56)
        )
        for sign in (1, -1):
            spins, energies = climb_batch(sign * matrix, initial)
            index = int(np.argmax(energies))
            value = int(energies[index])
            trials += batch_size
            if value > best_energy:
                best_energy = value
                best_spin = spins[index].copy()
                if sign < 0:
                    value = -value
            else:
                value = sign * value
            if best_energy >= threshold:
                actual = exact_energy(matrix, best_spin)
                if abs(actual) != best_energy:
                    raise AssertionError((actual, best_energy))
                return best_spin, actual, trials
    if best_spin is not None:
        actual = exact_energy(matrix, best_spin)
        if abs(actual) != best_energy:
            raise AssertionError((actual, best_energy))
    return None, None, trials


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source",
        type=Path,
        default=Path("computations/results/heuristic_m14_from_conference.json"),
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--threshold", type=int, default=210)
    parser.add_argument(
        "--diagonal-family",
        choices=("balanced", "all"),
        default="balanced",
    )
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--batches", type=int, default=4)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument(
        "--retry-holdouts",
        type=int,
        default=8,
        help="extra batches for first-pass holdouts",
    )
    parser.add_argument(
        "--toggle-micro-position",
        action="append",
        default=[],
        metavar="A,B",
        help="toggle position (A,B) in every ordered upper macro block",
    )
    args = parser.parse_args()

    source = json.loads(args.source.read_text())
    base = np.asarray(source["matrix"], dtype=np.int64)
    if base.shape != (14, 14):
        raise AssertionError("expected order-14 base")
    if not np.array_equal(base @ base, 13 * np.eye(14, dtype=np.int64)):
        raise AssertionError("conference identity failed")

    rng = np.random.default_rng(args.seed)
    toggles = tuple(
        tuple(int(value) for value in item.split(","))
        for item in args.toggle_micro_position
    )
    if any(len(item) != 2 for item in toggles):
        raise ValueError("micro position must have form A,B")
    started = time.time()
    successes = []
    holdouts: list[tuple[int, np.ndarray, int]] = []
    diagonals = candidate_diagonals(args.diagonal_family)
    for number, (mask, diagonal) in enumerate(diagonals, 1):
        matrix = lift(base, diagonal, toggles)
        spin, energy, trials = search_witness(
            matrix, rng, args.threshold, args.batches, args.batch_size
        )
        if spin is None:
            holdouts.append((mask, diagonal, trials))
        else:
            successes.append(
                {
                    "diagonal_positive_mask": f"{mask:04x}",
                    "energy": energy,
                    "spin_bits_little_endian": encode_spin(spin),
                    "trials": trials,
                }
            )
        if number % 250 == 0:
            print(
                f"first pass {number}/{len(diagonals)} "
                f"successes={len(successes)} holdouts={len(holdouts)}",
                flush=True,
            )

    unresolved = []
    for number, (mask, diagonal, first_trials) in enumerate(holdouts, 1):
        matrix = lift(base, diagonal, toggles)
        spin, energy, trials = search_witness(
            matrix,
            rng,
            args.threshold,
            args.retry_holdouts,
            args.batch_size,
        )
        if spin is None:
            unresolved.append(
                {
                    "diagonal_positive_mask": f"{mask:04x}",
                    "trials": first_trials + trials,
                }
            )
        else:
            successes.append(
                {
                    "diagonal_positive_mask": f"{mask:04x}",
                    "energy": energy,
                    "spin_bits_little_endian": encode_spin(spin),
                    "trials": first_trials + trials,
                }
            )
        if number % 100 == 0:
            print(
                f"retry {number}/{len(holdouts)} unresolved={len(unresolved)}",
                flush=True,
            )

    successes.sort(key=lambda row: row["diagonal_positive_mask"])
    # Independent exact verification of every saved compact witness.
    record_hash = hashlib.sha256()
    minimum = None
    for row in successes:
        mask = int(row["diagonal_positive_mask"], 16)
        diagonal = np.asarray(
            [1 if mask & (1 << i) else -1 for i in range(14)],
            dtype=np.int64,
        )
        spin = decode_spin(row["spin_bits_little_endian"])
        checked = exact_energy(lift(base, diagonal, toggles), spin)
        if checked != int(row["energy"]) or abs(checked) < args.threshold:
            raise AssertionError((row, checked))
        minimum = abs(checked) if minimum is None else min(minimum, abs(checked))
        record_hash.update(
            (
                row["diagonal_positive_mask"]
                + ":"
                + row["spin_bits_little_endian"]
                + ":"
                + str(row["energy"])
                + "\n"
            ).encode()
        )

    payload = {
        "schema": "quadratic-signing-diagonal-family-audit-v2",
        "classification": (
            "exhaustive enumeration of balanced diagonals; each successful "
            "lower bound is an exact Boolean witness; unresolved heuristic "
            "holdouts, if any, remain open"
        ),
        "source": str(args.source),
        "source_matrix_sha256": matrix_hash(base),
        "definition": (
            "S_D=A tensor H_4 + D tensor "
            "(H_4-diag(H_4))"
        ),
        "diagonal_family": args.diagonal_family,
        "cross_micro_position_toggles": [list(item) for item in toggles],
        "threshold": args.threshold,
        "strict_conference_threshold": 0.5 * 56**1.5,
        "seed": args.seed,
        "first_pass_batches": args.batches,
        "retry_batches": args.retry_holdouts,
        "batch_size": args.batch_size,
        "diagonal_count": len(diagonals),
        "certified_diagonal_count": len(successes),
        "unresolved_count": len(unresolved),
        "minimum_certified_absolute_energy": minimum,
        "canonical_record_sha256": record_hash.hexdigest(),
        "elapsed_seconds": time.time() - started,
        "records": successes,
        "unresolved": unresolved,
    }
    if args.diagonal_family == "balanced":
        payload["balanced_diagonal_count"] = len(diagonals)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(
        f"done certified={len(successes)}/{len(diagonals)} "
        f"unresolved={len(unresolved)} min={minimum} "
        f"elapsed={payload['elapsed_seconds']:.3f}s",
        flush=True,
    )
    print(f"wrote {args.output}")
    return 0 if not unresolved else 2


if __name__ == "__main__":
    raise SystemExit(main())
