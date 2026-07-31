#!/usr/bin/env python3
"""Independently verify the exact micro-position certificate-hitting result."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def hadamard4() -> np.ndarray:
    return np.asarray(
        [[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]],
        dtype=np.int64,
    )


def decode_spin(encoded: str) -> np.ndarray:
    bits = int(encoded, 16)
    return np.asarray(
        [1 if bits & (1 << i) else -1 for i in range(56)], dtype=np.int64
    ).reshape(14, 4)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--certificate",
        type=Path,
        default=Path("computations/results/phase2e_all_diagonal_family_audit.json"),
    )
    parser.add_argument(
        "--result",
        type=Path,
        default=Path("computations/results/phase2f_cross_fiber_micro.json"),
    )
    parser.add_argument(
        "--source",
        type=Path,
        default=Path("computations/results/heuristic_m14_from_conference.json"),
    )
    args = parser.parse_args()

    certificate = json.loads(args.certificate.read_text())
    result = json.loads(args.result.read_text())
    base = np.asarray(json.loads(args.source.read_text())["matrix"], dtype=np.int64)
    records = certificate["records"]
    if len(records) != 16384:
        raise AssertionError("incomplete all-diagonal certificate")
    energies = np.asarray([int(row["energy"]) for row in records], dtype=np.int64)
    h4 = hadamard4()
    contributions = np.empty((len(records), 16), dtype=np.int64)
    for row, record in enumerate(records):
        spin = decode_spin(str(record["spin_bits_little_endian"]))
        column = 0
        for a in range(4):
            for b in range(4):
                contributions[row, column] = sum(
                    int(base[i, j] * h4[a, b] * spin[i, a] * spin[j, b])
                    for i in range(14)
                    for j in range(i + 1, 14)
                )
                column += 1

    solved = result["results"]["micro"]
    if solved["status"] != "OPTIMAL" or solved["objective"] != 2:
        raise AssertionError("unexpected solver result")
    chosen = [int(column) for column in solved["chosen_columns"]]
    if chosen != [1, 4]:
        raise AssertionError(chosen)
    corrected = energies - 2 * np.sum(contributions[:, chosen], axis=1)
    if int(np.min(corrected)) != 122 or int(np.max(corrected)) != 176:
        raise AssertionError((int(np.min(corrected)), int(np.max(corrected))))
    if np.any(np.abs(corrected) > 208):
        raise AssertionError("chosen pair does not hit every saved witness")
    digest = hashlib.sha256(corrected.astype("<i8").tobytes()).hexdigest()
    if digest != solved["corrected_energy_sha256"]:
        raise AssertionError("corrected-energy hash mismatch")

    # An independent exhaustive lower bound needs only the empty correction
    # and the 16 singleton bundles.  Every one leaves an explicit violation,
    # while the chosen pair is feasible, proving the exact minimum is two.
    singleton_maxima = []
    candidates = [energies] + [
        energies - 2 * contributions[:, column] for column in range(16)
    ]
    for candidate in candidates:
        singleton_maxima.append(int(np.max(np.abs(candidate))))
    if min(singleton_maxima) <= 208:
        raise AssertionError("a correction with fewer than two bundles works")

    print(
        "verified exact micro-bundle optimum 2 by exhaustive size-0/1 lower "
        "bound and full chosen-pair feasibility"
    )
    print(
        f"chosen positions=[(0,1),(1,0)] corrected range=[122,176] hash={digest}"
    )
    print(f"empty/singleton maximum-|energy| values={singleton_maxima}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
