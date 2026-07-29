#!/usr/bin/env python3
"""Search symmetric circulant signings whose square has sparse support."""

from __future__ import annotations

import itertools
import math
import numpy as np


def symmetric_sequences(order: int):
    assert order % 2 == 0
    half = order // 2
    for bits in itertools.product((-1, 1), repeat=half):
        aa = np.empty(order, dtype=np.int8)
        aa[0] = 0
        for kk in range(1, half):
            aa[kk] = aa[order - kk] = bits[kk - 1]
        aa[half] = bits[half - 1]
        yield aa


def periodic_autocorrelation(aa: np.ndarray) -> np.ndarray:
    return np.array(
        [int(np.dot(aa.astype(int), np.roll(aa.astype(int), -dd)))
         for dd in range(len(aa))],
        dtype=int,
    )


def connected_support(order: int, corr: np.ndarray) -> bool:
    gg = order
    for dd in np.flatnonzero(corr[1:] != 0) + 1:
        gg = math.gcd(gg, int(dd))
    return gg == 1


def best_connected(order: int):
    best = None
    for aa in symmetric_sequences(order):
        corr = periodic_autocorrelation(aa)
        if not connected_support(order, corr):
            continue
        degree = int(np.count_nonzero(corr[1:]))
        if best is None or degree < best[0]:
            best = (degree, aa.copy(), corr.copy())
    return best


if __name__ == "__main__":
    for nn in range(4, 31, 2):
        degree, seq, corr = best_connected(nn)
        symbols = "".join("0" if x == 0 else "+" if x == 1 else "-" for x in seq)
        print(nn, degree, symbols, corr.tolist())
