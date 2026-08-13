#!/usr/bin/env python3
"""Exhaustively check the asymmetric-to-common-spin recoupling certificates.

All energies in this script use the doubled normalization z.T @ A @ z.
The script checks every saved exact minimizer orbit representative of orders
3 through 8, the saved exact order-10 representative, and symmetric
conference matrices of orders 6, 10, 14, and 18.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
RESULTS = ROOT / "computations" / "results"
ARCSIN_KAPPA = math.pi / 2 - 1


def projective_spins(n: int) -> np.ndarray:
    """All Boolean spins modulo global negation, with coordinate zero fixed."""
    masks = np.arange(1 << (n - 1), dtype=np.uint64)[:, None]
    bits = (
        (masks >> np.arange(n - 1, dtype=np.uint64)[None, :]) & 1
    ).astype(np.int8)
    return np.concatenate(
        (np.ones((len(masks), 1), dtype=np.int8), 1 - 2 * bits), axis=1
    )


def one_sided_caps(matrix: np.ndarray) -> tuple[int, int]:
    """Return (maximum energy, minus minimum energy)."""
    m = len(matrix)
    if m <= 1:
        return 0, 0
    spins = projective_spins(m).astype(np.int16)
    energies = np.einsum(
        "bi,ij,bj->b", spins, matrix, spins, optimize=True
    )
    return int(energies.max()), int(-energies.min())


def nuclear_floor(matrix: np.ndarray) -> float:
    """Common lower bound on each of the two one-sided Boolean caps."""
    m = len(matrix)
    if m <= 1:
        return 0.0
    nuclear = float(np.abs(np.linalg.eigvalsh(matrix)).sum())
    return max(0.0, nuclear / math.pi - (1 - 2 / math.pi) * m)


def projector_floor(matrix: np.ndarray, target_sign: int) -> float:
    """Sign-specific spectral-projector/Hermite lower certificate.

    target_sign=+1 certifies max_z z^T matrix z.
    target_sign=-1 certifies -min_z z^T matrix z.
    """
    m = len(matrix)
    if m <= 1:
        return 0.0
    eigenvalues = np.linalg.eigvalsh(target_sign * matrix)
    positive = eigenvalues[eigenvalues > 1e-9]
    rank = len(positive)
    if rank == 0:
        return 0.0
    mass = float(positive.sum())
    theta = min(1.0, mass / (2 * ARCSIN_KAPPA * rank))
    return max(
        0.0,
        (2 / math.pi)
        * (theta * mass - ARCSIN_KAPPA * theta * theta * rank),
    )


def aligned_cap(
    positive_cap: int, negative_cap: int, target_sign: int
) -> int:
    return positive_cap if target_sign > 0 else negative_cap


def bilinear_maximizers(
    matrix: np.ndarray,
) -> tuple[int, list[tuple[np.ndarray, np.ndarray]]]:
    """Enumerate all projective x and all y attaining max x^T A y."""
    n = len(matrix)
    spins = projective_spins(n).astype(np.int16)
    fields = spins @ matrix
    row_values = np.abs(fields).sum(axis=1)
    bilinear_cap = int(row_values.max())
    pairs: list[tuple[np.ndarray, np.ndarray]] = []
    for index in np.flatnonzero(row_values == bilinear_cap):
        x = spins[index]
        field = fields[index]
        zero_coordinates = np.flatnonzero(field == 0)
        base_y = np.where(field >= 0, 1, -1).astype(np.int16)
        for mask in range(1 << len(zero_coordinates)):
            y = base_y.copy()
            for bit, coordinate in enumerate(zero_coordinates):
                y[coordinate] = -1 if (mask >> bit) & 1 else 1
            assert int(x @ matrix @ y) == bilinear_cap
            pairs.append((x.copy(), y))
    return bilinear_cap, pairs


def audit_matrix(matrix: np.ndarray) -> dict[str, object]:
    n = len(matrix)
    spins = projective_spins(n).astype(np.int16)
    energies = np.einsum(
        "bi,ij,bj->b", spins, matrix, spins, optimize=True
    )
    quadratic_cap = int(np.abs(energies).max())
    bilinear_cap, pairs = bilinear_maximizers(matrix)

    best_exact = -math.inf
    best_projector = -math.inf
    best_nuclear = -math.inf
    best_data: tuple[int, int, int, int] | None = None

    for x, y in pairs:
        agreement = np.flatnonzero(x == y)
        disagreement = np.flatnonzero(x != y)
        block_i = matrix[np.ix_(agreement, agreement)]
        block_j = matrix[np.ix_(disagreement, disagreement)]
        p = x[agreement]
        q = x[disagreement]
        energy_i = int(p @ block_i @ p) if len(agreement) else 0
        energy_j = int(q @ block_j @ q) if len(disagreement) else 0
        assert energy_i - energy_j == bilinear_cap

        if energy_i * energy_j >= 0:
            exact = projector = nuclear = float(abs(energy_i - energy_j))
        else:
            sign_i = 1 if energy_i > 0 else -1
            sign_j = 1 if energy_j > 0 else -1
            pos_i, neg_i = one_sided_caps(block_i)
            pos_j, neg_j = one_sided_caps(block_j)

            exact = float(
                max(
                    abs(energy_i)
                    + aligned_cap(pos_j, neg_j, sign_i),
                    abs(energy_j)
                    + aligned_cap(pos_i, neg_i, sign_j),
                )
            )
            projector = max(
                abs(energy_i) + projector_floor(block_j, sign_i),
                abs(energy_j) + projector_floor(block_i, sign_j),
            )
            nuclear = max(
                abs(energy_i) + nuclear_floor(block_j),
                abs(energy_j) + nuclear_floor(block_i),
            )

        assert exact <= quadratic_cap + 1e-8
        assert projector <= exact + 1e-8
        assert nuclear <= exact + 1e-8
        if exact > best_exact:
            best_exact = exact
            best_data = (
                len(agreement),
                len(disagreement),
                energy_i,
                energy_j,
            )
        if projector > best_projector:
            best_projector = projector
        if nuclear > best_nuclear:
            best_nuclear = nuclear

    return {
        "n": n,
        "quadratic_cap": quadratic_cap,
        "bilinear_cap": bilinear_cap,
        "polarization": bilinear_cap / 2,
        "exact_recoupling": best_exact,
        "projector_recoupling": best_projector,
        "nuclear_recoupling": best_nuclear,
        "best_pair": best_data,
        "bilinear_maximizer_count": len(pairs),
    }


def load_cases() -> list[tuple[str, np.ndarray]]:
    cases: list[tuple[str, np.ndarray]] = []
    for n in range(3, 9):
        payload = json.loads(
            (RESULTS / f"m{n}_minimizer_orbits.json").read_text()
        )
        for row in payload["classes"]:
            cases.append(
                (
                    f"min-n{n}-orbit{row['class']}",
                    np.asarray(row["representative_matrix"], dtype=np.int16),
                )
            )

    exact_ten = json.loads((RESULTS / "exact_m10.json").read_text())
    cases.append(
        ("exact-n10", np.asarray(exact_ten["matrix"], dtype=np.int16))
    )

    conference_sources = (
        ("conference-n6", "conference_double_p5.json"),
        ("conference-n10", "conference_order10_gf9.json"),
        ("conference-n14", "conference_double_p13.json"),
        ("conference-n18", "conference_double_p17.json"),
    )
    for label, filename in conference_sources:
        payload = json.loads((RESULTS / filename).read_text())
        cases.append(
            (label, np.asarray(payload["conference_matrix"], dtype=np.int16))
        )
    return cases


def main() -> int:
    print(
        "| case | n | Q | B | B/2 | exact face | projector | nuclear "
        "| best (|I|,|J|,P,R) |"
    )
    print("|---|---:|---:|---:|---:|---:|---:|---:|---|")
    for label, matrix in load_cases():
        row = audit_matrix(matrix)
        print(
            f"| {label} | {row['n']} | {row['quadratic_cap']} "
            f"| {row['bilinear_cap']} | {row['polarization']:.3f} "
            f"| {row['exact_recoupling']:.3f} "
            f"| {row['projector_recoupling']:.3f} "
            f"| {row['nuclear_recoupling']:.3f} "
            f"| {row['best_pair']} |"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
