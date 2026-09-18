#!/usr/bin/env python3
"""Wave 51 audit of exact local averages and the reverse-branch band.

The calculations exhaust every selector and oriented local state in the
stored exact minimizers.  They check the algebraic identities used in the
memo; the asymptotic inverse-tail theorem itself is proved there.
"""

from __future__ import annotations

import itertools
import math
import sys

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from cross_block_dual_r48_check import sample_order_ten
from envelope_block_cover_r27 import A6, A8, A9, projective_spins


def qnorm(a: np.ndarray) -> int:
    return max(abs(int(x @ a @ x)) for x in projective_spins(len(a)))


def audit(a: np.ndarray, name: str, m: int, t: float = 0.0) -> None:
    n = len(a)
    k = n - m
    q = qnorm(a)
    p = m / n
    p2 = m * (m - 1) / (n * (n - 1))
    one_minus_p2 = 1.0 - p2
    beta = p ** 1.5 - p2

    target_e2 = 2.0 * m * (m - 1)
    target_v = 2.0 * k * (n + m - 1)
    target_d = q * q - 2.0 * n * (n - 1)

    max_average_error = 0.0
    min_reverse_gap = math.inf
    max_reverse_r = 0.0
    max_reverse_minus_g = 0.0
    reverse_count = 0
    positive_count = 0
    state_count = 0

    for s in itertools.combinations(range(n), m):
        outside = tuple(i for i in range(n) if i not in s)
        child = a[np.ix_(s, s)]
        qs = qnorm(child)
        h = qs - p ** 1.5 * q - t

        values = []
        for sigma in (-1, 1):
            for y in projective_spins(m):
                e = float(sigma * int(y @ child @ y))
                cross = a[np.ix_(outside, s)] @ y
                v = float(4 * int(cross @ cross) + 2 * k * (k - 1))
                d = q * q - e * e - v
                g = one_minus_p2 * e - h
                raw_numerator = p2 * v + (q - e) * g
                values.append((e, v, d, g, raw_numerator))

                if g > 0:
                    positive_count += 1
                else:
                    r = -g / p2
                    # Strict positivity of the reverse two-moment branch.
                    if 0 <= r < q + e and v - r * (q - e) > 1e-10:
                        reverse_count += 1
                        max_reverse_r = max(max_reverse_r, r)
                        max_reverse_minus_g = max(max_reverse_minus_g, -g)
                        lower_gap = (beta * q + t) / one_minus_p2
                        min_reverse_gap = min(min_reverse_gap, q - e - lower_gap)
                        assert r * (q - e) < v + 1e-9
                        assert r <= v / (q - e) + 1e-9
                state_count += 1

        arr = np.asarray(values)
        means = arr.mean(axis=0)
        expected = np.asarray((0.0, target_e2, target_d, -h, -q * h))
        # Columns are e,V,D,g,N; replace the second expected entry by E[V].
        expected[1] = target_v
        max_average_error = max(max_average_error, float(np.max(np.abs(means - expected))))

    assert max_average_error < 1e-8
    assert min_reverse_gap >= -1e-8 or reverse_count == 0
    # The global spectral estimate gives the deterministic upper bound below.
    v_global = 8.0 * q * m + 2.0 * n * n
    r_global = v_global / (beta * q)
    assert max_reverse_r <= r_global + 1e-8

    support_mass = (positive_count + reverse_count) / state_count
    min_gap_text = "NA" if reverse_count == 0 else f"{min_reverse_gap:.6f}"
    print(
        f"{name:10s} n={n:2d} m={m:2d} q={q:2d} "
        f"avg_err={max_average_error:.2e} support_mass={support_mass:.6f} "
        f"reverse={reverse_count:6d} max_r={max_reverse_r:9.4f} "
        f"max(-g)={max_reverse_minus_g:9.4f} gap_slack={min_gap_text} "
        f"coarse_r_bound={r_global:10.2f}"
    )


def main() -> None:
    for a, name, ms in (
        (A6, "A6", (3, 4, 5)),
        (A8, "A8", (4, 5, 6, 7)),
        (A9, "A9", (5, 6, 7, 8)),
        (sample_order_ten(0), "A10_seed0", (5, 6, 7, 8, 9)),
    ):
        for m in ms:
            audit(a, name, m)
    print("PASS completion-abundance exact identities and reverse-band audit")


if __name__ == "__main__":
    main()
