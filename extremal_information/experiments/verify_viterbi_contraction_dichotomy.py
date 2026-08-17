#!/usr/bin/env python3
"""Exact 2x2 checks for the max-plus contraction dichotomy."""

from __future__ import annotations

from fractions import Fraction
from itertools import product
import json


def update(a, K):
    return tuple(max(a[i] + K[i][j] for i in range(2)) for j in range(2))


def span(z):
    return max(z) - min(z)


def subtract(a, b):
    return tuple(x - y for x, y in zip(a, b))


def rank_one(K):
    return K[0][0] + K[1][1] == K[0][1] + K[1][0]


def witness_ratio_one(K):
    d0 = K[0][0] - K[1][0]
    d1 = K[0][1] - K[1][1]
    if d0 == d1:
        return False
    r = Fraction(d0 + d1, 2)
    a = (Fraction(0), r)
    gap = abs(d0 - d1) / 2
    delta = Fraction(1, 2) * gap
    b = (delta, r)
    input_span = span(subtract(a, b))
    output_span = span(subtract(update(a, K), update(b, K)))
    return input_span > 0 and output_span == input_span


def verify():
    exact_resets = 0
    exact_ratio_one = 0
    for entries in product(range(-2, 3), repeat=4):
        K = (entries[:2], entries[2:])
        if rank_one(K):
            outputs = set()
            for a in product(range(-3, 4), repeat=2):
                out = update(a, K)
                outputs.add(tuple(value - max(out) for value in out))
            assert len(outputs) == 1
            exact_resets += 1
        else:
            assert witness_ratio_one(K)
            exact_ratio_one += 1
    return {"rank_one_reset_matrices": exact_resets,
            "nonrank_one_ratio_one_witnesses": exact_ratio_one,
            "matrix_entry_range": [-2, 2],
            "arithmetic": "exact rational"}


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2, sort_keys=True))
