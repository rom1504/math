#!/usr/bin/env python3
"""Exact finite checks for the regular-Hadamard local-stability obstruction."""

from __future__ import annotations

import argparse
import json

import numpy as np


def regular_hadamard(k: int) -> tuple[np.ndarray, np.ndarray]:
    s = 1 << k
    m = s * s
    values = np.arange(m, dtype=np.int64)
    mask = s - 1
    a = values & mask
    b = values >> k
    parity = np.array([bin(i).count("1") & 1 for i in range(s)], dtype=np.int8)
    bent_sign = 1 - 2 * parity[a & b]
    matrix = bent_sign[np.bitwise_xor(values[:, None], values[None, :])]
    in_subspace = (
        (((b >> 1) & 1) == ((a >> 1) & 1))
        & ((b & 1) == ((a & 1) ^ ((a >> 1) & 1)))
    )
    return matrix.astype(np.int32), in_subspace


def best_improvement(c_matrix: np.ndarray) -> tuple[np.ndarray, int, int]:
    """Maximize -r^T C r from r=1, breaking gain ties by first index."""
    spin = np.ones(len(c_matrix), dtype=np.int32)
    field = c_matrix @ spin
    flips = 0
    while True:
        gains = 4 * spin * field
        coordinate = int(np.argmax(gains))
        if int(gains[coordinate]) <= 0:
            return spin, int(-(spin @ c_matrix @ spin)), flips
        old = int(spin[coordinate])
        spin[coordinate] = -old
        field -= 2 * old * c_matrix[:, coordinate]
        flips += 1


def check(k: int) -> dict[str, object]:
    hadamard, in_subspace = regular_hadamard(k)
    m = len(hadamard)
    s = 1 << k
    identity = np.eye(m, dtype=np.int32)
    ones = np.ones(m, dtype=np.int32)
    c_matrix = hadamard - identity
    stable_spin = np.where(in_subspace, -1, 1).astype(np.int32)

    assert np.array_equal(hadamard @ hadamard, m * identity)
    assert np.array_equal(hadamard @ ones, s * ones)
    assert int(in_subspace.sum()) == m // 4
    assert np.array_equal(
        hadamard @ in_subspace.astype(np.int32),
        (s // 2) * stable_spin,
    )
    products = stable_spin * (c_matrix @ stable_spin)
    assert np.all(products[~in_subspace] == -1)
    assert np.all(products[in_subspace] == -(2 * s + 1))
    explicit_energy = int(-(stable_spin @ c_matrix @ stable_spin))
    assert explicit_energy == m * (s // 2 + 1)

    terminal_spin, terminal_energy, flips = best_improvement(c_matrix)
    return {
        "k": k,
        "m": m,
        "s": s,
        "explicit_stable_energy": explicit_energy,
        "greedy_terminal_energy": terminal_energy,
        "greedy_flip_count": flips,
        "greedy_reaches_explicit_spin": bool(np.array_equal(terminal_spin, stable_spin)),
        "normalized_possible_joint_defect": (m * (s / 2 - 2))
        / ((2 * m) ** 1.5),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-k", type=int, default=4)
    args = parser.parse_args()
    if args.max_k < 2:
        raise ValueError("max-k must be at least 2")
    print(json.dumps([check(k) for k in range(2, args.max_k + 1)], indent=2))


if __name__ == "__main__":
    main()
