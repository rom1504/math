#!/usr/bin/env python3
"""Independently verify the balanced-diagonal family witness certificate."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from collections import Counter
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


def expected_masks(family: str) -> set[int]:
    if family == "balanced":
        return {
            sum(1 << vertex for vertex in positive)
            for positive in itertools.combinations(range(14), 7)
        }
    if family == "all":
        return set(range(1 << 14))
    raise AssertionError(f"unknown diagonal family: {family}")


def decode_spin(encoded: str) -> np.ndarray:
    if len(encoded) != 14:
        raise AssertionError("a 56-bit spin must use exactly 14 hex digits")
    bits = int(encoded, 16)
    return np.asarray(
        [1 if bits & (1 << i) else -1 for i in range(56)], dtype=np.int64
    )


def lift(
    base: np.ndarray,
    mask: int,
    require_balanced: bool,
    cross_micro_position_toggles: tuple[tuple[int, int], ...] = (),
) -> np.ndarray:
    diagonal = np.asarray(
        [1 if mask & (1 << i) else -1 for i in range(14)], dtype=np.int64
    )
    if require_balanced and int(diagonal.sum()) != 0:
        raise AssertionError("diagonal is not balanced")
    hadamard = sylvester4()
    off_hadamard = hadamard - np.diag(np.diag(hadamard))
    matrix = np.kron(base, hadamard) + np.kron(
        np.diag(diagonal), off_hadamard
    )
    for a, b in cross_micro_position_toggles:
        if not (0 <= a < 4 and 0 <= b < 4):
            raise AssertionError((a, b))
        for i in range(14):
            for j in range(i + 1, 14):
                u = 4 * i + a
                v = 4 * j + b
                matrix[u, v] *= -1
                matrix[v, u] *= -1
    if not np.array_equal(matrix, matrix.T) or np.any(np.diag(matrix)):
        raise AssertionError("invalid lifted signing")
    if not np.all(np.abs(matrix[~np.eye(56, dtype=bool)]) == 1):
        raise AssertionError("invalid lifted off-diagonal entry")
    return matrix


def exact_energy(matrix: np.ndarray, spin: np.ndarray) -> int:
    quadratic = int(spin @ matrix @ spin)
    if quadratic % 2:
        raise AssertionError("odd quadratic value")
    return quadratic // 2


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "certificate",
        type=Path,
        nargs="?",
        default=Path(
            "computations/results/phase2d_balanced_diagonal_family_audit.json"
        ),
    )
    parser.add_argument(
        "--source",
        type=Path,
        default=Path("computations/results/heuristic_m14_from_conference.json"),
    )
    args = parser.parse_args()

    certificate = json.loads(args.certificate.read_text())
    source = json.loads(args.source.read_text())
    base = np.asarray(source["matrix"], dtype=np.int64)
    if base.shape != (14, 14):
        raise AssertionError("expected order-14 base")
    if not np.array_equal(base, base.T) or np.any(np.diag(base)):
        raise AssertionError("invalid base signing")
    if not np.array_equal(base @ base, 13 * np.eye(14, dtype=np.int64)):
        raise AssertionError("base conference identity failed")
    if matrix_hash(base) != certificate["source_matrix_sha256"]:
        raise AssertionError("source matrix hash mismatch")
    if certificate["unresolved"] or certificate["unresolved_count"] != 0:
        raise AssertionError("certificate still has unresolved diagonals")

    records = certificate["records"]
    masks = [int(record["diagonal_positive_mask"], 16) for record in records]
    family = certificate.get("diagonal_family", "balanced")
    toggles = tuple(
        tuple(int(value) for value in item)
        for item in certificate.get("cross_micro_position_toggles", [])
    )
    expected = expected_masks(family)
    expected_count = len(expected)
    if (
        len(records) != expected_count
        or set(masks) != expected
        or len(set(masks)) != expected_count
    ):
        raise AssertionError("records do not enumerate the requested diagonals once")
    recorded_count = certificate.get(
        "diagonal_count", certificate.get("balanced_diagonal_count")
    )
    if recorded_count != expected_count:
        raise AssertionError("incorrect diagonal count")
    if certificate["certified_diagonal_count"] != expected_count:
        raise AssertionError("incorrect certified count")

    threshold = int(certificate["threshold"])
    if threshold != 210:
        raise AssertionError("unexpected threshold")
    # Exact strict comparison threshold / 56^(3/2) > 1/2.
    if 4 * threshold * threshold <= 56**3:
        raise AssertionError("threshold does not strictly exceed constant 1/2")

    canonical_hash = hashlib.sha256()
    energies = Counter()
    minimum = None
    for record in sorted(records, key=lambda row: row["diagonal_positive_mask"]):
        mask_text = record["diagonal_positive_mask"]
        spin_text = record["spin_bits_little_endian"]
        claimed = int(record["energy"])
        matrix = lift(
            base,
            int(mask_text, 16),
            family == "balanced",
            toggles,
        )
        checked = exact_energy(matrix, decode_spin(spin_text))
        if checked != claimed or abs(checked) < threshold:
            raise AssertionError((mask_text, claimed, checked))
        energies[checked] += 1
        minimum = abs(checked) if minimum is None else min(minimum, abs(checked))
        canonical_hash.update(
            (mask_text + ":" + spin_text + ":" + str(claimed) + "\n").encode()
        )

    if minimum != certificate["minimum_certified_absolute_energy"]:
        raise AssertionError("minimum energy mismatch")
    if canonical_hash.hexdigest() != certificate["canonical_record_sha256"]:
        raise AssertionError("canonical record hash mismatch")

    # For the uncorrected common-Sylvester family, verify the tensor channel:
    # v is Boolean and H_4 v = 2v.  A nonempty micro-position correction is
    # only a finite order-56 certificate and does not inherit this claim.
    hadamard = sylvester4()
    micro = np.asarray([-1, -1, -1, 1], dtype=np.int64)
    if int(np.trace(hadamard)) != 0:
        raise AssertionError("order-four Hadamard diagonal is not balanced")
    if not np.array_equal(hadamard @ micro, 2 * micro):
        raise AssertionError("micro Boolean eigenvector identity failed")

    print(
        f"verified all {expected_count} {family} diagonals exactly; "
        f"minimum |energy|={minimum}; energy distribution={dict(sorted(energies.items()))}"
    )
    if not toggles:
        print(
            "every fixed-D Sylvester descendant has normalized cap at least "
            f"{threshold / 56**1.5:.15f} > 1/2"
        )
    else:
        print(
            "micro-position toggles are present: verified finite order-56 "
            "witnesses only; no tensor-persistence claim"
        )
    print(f"canonical records sha256={canonical_hash.hexdigest()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
