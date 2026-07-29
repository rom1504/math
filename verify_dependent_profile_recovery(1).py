#!/usr/bin/env python3
"""Exhaustively verify the dependent 4-fibre compressed-lift witness."""

import json
from itertools import product
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent


def main():
    data = json.loads(
        (HERE / "dependent_profile_recovery_witness.json").read_text()
    )
    matrix = np.asarray(data["matrix"], dtype=np.int64)
    n = int(data["seed_order"])
    s = int(data["fibre_size"])
    order = n * s

    assert matrix.shape == (order, order)
    assert np.array_equal(matrix, matrix.T)
    assert np.all(np.diag(matrix) == 0)
    assert set(matrix[~np.eye(order, dtype=bool)]) == {-1, 1}

    block_sums = [
        int(matrix[i * s : (i + 1) * s, j * s : (j + 1) * s].sum())
        for i, j in data["seed_edge_order"]
    ]
    assert block_sums == data["required_cross_block_sums"]

    # Gauge x_0=1 because x and -x have identical quadratic energy.
    best_abs = -1
    best_energy = None
    best_spin = None
    histogram = {}
    for tail in product((-1, 1), repeat=order - 1):
        spin = np.asarray((1,) + tail, dtype=np.int64)
        energy = int(spin @ matrix @ spin)
        histogram[energy] = histogram.get(energy, 0) + 1
        if abs(energy) > best_abs:
            best_abs = abs(energy)
            best_energy = energy
            best_spin = tuple(map(int, spin))

    assert best_abs == data["boolean_quadratic_norm"] == 40
    assert best_abs / order ** 1.5 == data["doubled_normalized_value"]
    assert best_abs / (2 * order ** 1.5) == data[
        "original_half_energy_normalized_value"
    ]
    assert data["uncontracted_target"] == s ** 1.5 * data[
        "seed_boolean_quadratic_norm"
    ]

    # Verify the exact chiral signed antiautomorphism.
    permutation = data["signed_antiautomorphism"]["permutation_image"]
    signs = data["signed_antiautomorphism"]["column_signs"]
    signed_permutation = np.zeros((order, order), dtype=np.int64)
    for source, (target, sign) in enumerate(zip(permutation, signs)):
        signed_permutation[target, source] = sign
    assert np.array_equal(
        signed_permutation @ signed_permutation,
        -np.eye(order, dtype=np.int64),
    )
    assert np.array_equal(
        signed_permutation.T @ matrix @ signed_permutation,
        -matrix,
    )
    assert np.array_equal(
        signed_permutation @ matrix,
        -matrix @ signed_permutation,
    )

    print("valid zero-diagonal signing of order", order)
    print("cross-block sums:", block_sums)
    print("exhaustive Q:", best_abs, "witness energy:", best_energy)
    print("one maximizing spin:", best_spin)
    print("energy support:", sorted(histogram))
    print(
        "normalized doubled / original:",
        best_abs / order ** 1.5,
        best_abs / (2 * order ** 1.5),
    )
    print("signed antiautomorphism S^2=-I and SB=-BS: verified")


if __name__ == "__main__":
    main()
