#!/usr/bin/env python3
"""Block-palindrome circulant family with finite-range large correlations."""

from __future__ import annotations

import numpy as np


def block_palindrome(blocks: int, block_length: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    eps = rng.choice((-1, 1), size=blocks)
    half = np.repeat(eps, block_length)
    order = 2 * len(half)
    seq = np.empty(order, dtype=np.int8)
    seq[0] = 0
    seq[1 : len(half) + 1] = half
    for kk in range(len(half) + 1, order):
        seq[kk] = seq[order - kk]
    return seq


def periodic_autocorrelation(seq: np.ndarray) -> np.ndarray:
    vec = seq.astype(np.int64)
    return np.array(
        [np.dot(vec, np.roll(vec, -dd)) for dd in range(len(vec))],
        dtype=np.int64,
    )


if __name__ == "__main__":
    for blocks in (100, 500, 2000):
        aa = block_palindrome(blocks, block_length=8, seed=blocks)
        rr = periodic_autocorrelation(aa) / (len(aa) - 1)
        print(
            "n", len(aa),
            "q[1:10]", np.round(rr[1:10], 4).tolist(),
            "mid-max", float(np.max(np.abs(rr[8:-8]))),
        )
