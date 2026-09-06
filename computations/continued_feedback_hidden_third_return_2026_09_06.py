#!/usr/bin/env python3
"""Actual signing stress test: flat B² off-diagonal need not flatten B³.

Randomness chooses explicit sign gauges only. Operator bound is deterministic.
The fixed three-block path identity is checked in integer arithmetic.
"""

import argparse
import json

import numpy as np

from continued_feedback_variance_falsifier_2026_09_06 import sylvester


def construct(m: int, seed: int):
    h = sylvester(m)
    root = int(np.sqrt(m))
    assert root * root == m
    bits = (m.bit_length() - 1)
    assert bits % 2 == 0
    bent = np.array([
        -1 if sum(((x >> j) & 1) * ((x >> (j + 1)) & 1)
                  for j in range(0, bits, 2)) % 2 else 1
        for x in range(m)
    ], dtype=np.int64)
    third_numerator = h @ (bent[:, None] * h)
    assert np.all(third_numerator % root == 0)
    third = third_numerator // root
    assert np.all(np.abs(third) == 1)
    middle = bent[:, None] * h
    assert np.array_equal(h @ middle @ third, m * root * np.eye(m, dtype=np.int64))
    rng = np.random.default_rng(seed)
    blocks = [[None] * 4 for _ in range(4)]
    for i in range(4):
        gauge = rng.choice((-1, 1), size=m)
        blocks[i][i] = gauge[:, None] * h * gauge[None, :]
    fixed = {(0, 1): h, (1, 2): middle, (2, 3): third}
    for i in range(4):
        for j in range(i + 1, 4):
            if (i, j) in fixed:
                block = fixed[i, j]
            else:
                left = rng.choice((-1, 1), size=m)
                right = rng.choice((-1, 1), size=m)
                block = left[:, None] * h * right[None, :]
            blocks[i][j] = block
            blocks[j][i] = block.T
    a = np.block(blocks)
    np.fill_diagonal(a, 0)
    assert np.array_equal(a, a.T)
    assert np.all(np.abs(a[np.triu_indices(4 * m, 1)]) == 1)
    return a


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sizes", nargs="+", type=int, default=[4, 16, 64, 256])
    parser.add_argument("--seed", type=int, default=8092026)
    args = parser.parse_args()
    for m in args.sizes:
        a = construct(m, args.seed + m)
        n = len(a)
        b = a.astype(float) / np.sqrt(n - 1)
        q = b @ b
        b3 = q @ b
        off = q - np.diag(q.diagonal())
        return_block = b3[:m, 3 * m :]
        diagonal = np.diag(return_block)
        record = {
            "status": "finite numerical stress test; fixed path identity exact",
            "m": m,
            "n": n,
            "seed": args.seed + m,
            "deterministic_operator_bound": float((4 * np.sqrt(m) + 1) / np.sqrt(4 * m - 1)),
            "max_offdiag_B2": float(np.max(np.abs(off))),
            "mean_target_B3_diagonal": float(np.mean(diagonal)),
            "max_target_error_from_one_eighth": float(np.max(np.abs(diagonal - 0.125))),
            "normalized_cubic_trace": float(np.trace(b3) / n),
        }
        print(json.dumps(record), flush=True)


if __name__ == "__main__":
    main()
