#!/usr/bin/env python3
"""Exact Walsh32 bilinear norm160 and regular-H4-square/H2 half-cap80.

Uses only integer enumeration of32768 Walsh16 inputs and one840×448
exceptional profile comparison. No optimization solver or floating point.
"""

from __future__ import annotations

import json

import numpy as np


def walsh(order: int) -> np.ndarray:
    result = np.ones((1, 1), dtype=np.int64)
    while len(result) < order:
        result = np.block([[result, result], [result, -result]])
    assert len(result) == order
    return result


def unrestricted_matching_bound(a: np.ndarray, b: np.ndarray) -> int:
    # Every coordinate chooses either its a or its b term in max(a_i,b_i).
    # If k terms choose a, their sum is at most the k largest a entries;
    # similarly the remaining16-k b terms. This does not assume realizability.
    a_prefix = np.concatenate(([0], np.cumsum(a[::-1])))
    b_prefix = np.concatenate(([0], np.cumsum(b[::-1])))
    return 2 * max(int(a_prefix[k] + b_prefix[16 - k]) for k in range(17))


def main() -> None:
    w16 = walsh(16)
    codes = np.arange(1 << 15, dtype=np.uint64)
    spins = 1 - 2 * ((codes[:, None] >> np.arange(16, dtype=np.uint64)) & 1).astype(np.int64)
    absolute = np.abs(spins @ w16)
    types, labels, counts = np.unique(
        np.sort(absolute, axis=1), axis=0, return_inverse=True, return_counts=True
    )
    assert len(types) == 8
    exceptional = []
    pair_upper = 0
    for i, a in enumerate(types):
        for j, b in enumerate(types):
            bound = unrestricted_matching_bound(a, b)
            if bound <= 160:
                pair_upper = max(pair_upper, bound)
                continue
            first = np.unique(absolute[labels == i], axis=0)
            second = np.unique(absolute[labels == j], axis=0)
            values = 2 * np.maximum(first[:, None, :], second[None, :, :]).sum(axis=2)
            exact_max = int(values.max())
            assert exact_max <= 160
            pair_upper = max(pair_upper, exact_max)
            levels, multiplicities = np.unique(values, return_counts=True)
            exceptional.append({
                "types": [i, j],
                "unrestricted_upper": bound,
                "profile_counts": [len(first), len(second)],
                "exact_pair_max": exact_max,
                "value_counts": [[int(x), int(y)] for x, y in zip(levels, multiplicities)],
            })
    assert pair_upper == 160

    # The block Walsh recurrence proves ||W32 v||1≤160 for everyBoolean v.
    # Independent row/column signs and row permutations preserve this norm.
    h4 = np.ones((4, 4), dtype=np.int64) - 2 * np.eye(4, dtype=np.int64)
    dr = np.diag((-1, 1, 1, 1))
    dc = np.diag((1, -1, -1, -1))
    assert np.array_equal(dr @ h4 @ dc, walsh(4)[[0, 2, 1, 3]])
    seed = np.array([[1, 1], [1, -1]], dtype=np.int64)
    matrix = np.kron(h4, np.kron(h4, seed))
    witness = np.array([
        1, 1, 1, 1, 1, -1, 1, 1,
        -1, -1, -1, -1, -1, 1, -1, -1,
        1, 1, 1, 1, 1, -1, 1, 1,
        -1, -1, -1, -1, -1, 1, -1, -1,
    ], dtype=np.int64)
    doubled_energy = int(witness @ matrix @ witness)
    assert doubled_energy == -160
    print(json.dumps({
        "classification": "exact exhaustive integer certificate",
        "projective_walsh16_inputs": len(codes),
        "absolute_spectrum_types": types.tolist(),
        "type_counts": counts.tolist(),
        "exceptional_type_pairs": exceptional,
        "walsh32_bilinear_norm": 160,
        "regular_h4_square_h2_absolute_half_cap": 80,
        "quadratic_witness": witness.tolist(),
        "normalized_regularized_seed_witness": "80 / 16^(3/2) = 5/4",
        "scope": "finite outer order16 only; not an upper bound forR",
    }, indent=2))


if __name__ == "__main__":
    main()
