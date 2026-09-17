#!/usr/bin/env python3
"""Check two exact one-sided lower bounds for the completion CDF."""

from __future__ import annotations

import itertools

import numpy as np

from cross_block_dual_r48_check import sample_order_ten
from envelope_block_cover_r27 import A6, A8, A9, projective_spins


def qnorm(a: np.ndarray) -> int:
    return max(abs(int(x @ a @ x)) for x in projective_spins(len(a)))


def audit(name: str, a_mat: np.ndarray, m: int) -> None:
    n = len(a_mat)
    q = qnorm(a_mat)
    p = m / n
    p2 = m * (m - 1) / (n * (n - 1))
    p32q = p ** 1.5 * q
    actual_sum = cantelli_sum = reverse_sum = combined_sum = 0.0
    states = 0
    positive_states = reverse_active_states = 0
    for s in itertools.combinations(range(n), m):
        s = tuple(s)
        tset = tuple(i for i in range(n) if i not in s)
        child = a_mat[np.ix_(s, s)]
        qs = qnorm(child)
        # Gauge-anchor the first selected coordinate.  Together with sigma,
        # these are the 2^m oriented local states.
        for tail in itertools.product((-1, 1), repeat=m - 1):
            y = np.asarray((1,) + tail, dtype=np.int64)
            for sigma in (-1, 1):
                e = sigma * int(y @ child @ y)
                g = (1 - p2) * e - (qs - p32q)
                zvals = []
                for outside in itertools.product((-1, 1), repeat=n - m):
                    x = np.empty(n, dtype=np.int64)
                    x[list(s)] = y
                    x[list(tset)] = outside
                    zvals.append(sigma * int(x @ a_mat @ x) - e)
                zvals = np.asarray(zvals, dtype=float)
                assert abs(float(np.mean(zvals))) < 1e-12
                v = float(np.mean(zvals**2))
                actual = float(np.mean(p2 * zvals <= g + 1e-12))
                cantelli = reverse = 0.0
                if g > 0:
                    positive_states += 1
                    cantelli = g * g / (g * g + p2 * p2 * v)
                elif g < 0:
                    r = -g / p2
                    aa = q + e
                    b = q - e
                    if 0 < r < aa and v > r * b:
                        reverse_active_states += 1
                        reverse = (v - r * b) / ((aa - r) * (aa + b))
                assert cantelli <= actual + 1e-9
                assert reverse <= actual + 1e-9
                actual_sum += actual
                cantelli_sum += cantelli
                reverse_sum += reverse
                combined_sum += cantelli + reverse
                states += 1
    print(
        name,
        "m", m,
        "Z0", round(actual_sum / states, 9),
        "Cantelli", round(cantelli_sum / states, 9),
        "reverse", round(reverse_sum / states, 9),
        "combined", round(combined_sum / states, 9),
        "capture", round(combined_sum / actual_sum, 6),
        "gpos_frac", round(positive_states / states, 6),
        "reverse_active_frac", round(reverse_active_states / states, 6),
    )


def main() -> None:
    for name, a, ms in [
        ("A6", A6, (3, 5)),
        ("A8", A8, (4, 7)),
        ("A9", A9, (5, 8)),
        ("A10s0", sample_order_ten(0), (5, 9)),
    ]:
        for m in ms:
            audit(name, a, m)


if __name__ == "__main__":
    main()
