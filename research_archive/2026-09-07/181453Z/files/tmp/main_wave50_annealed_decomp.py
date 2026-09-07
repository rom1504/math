#!/usr/bin/env python3
"""Independent Wave 50 diagnostic for the two sources of annealed mass."""

from __future__ import annotations

import itertools
from collections import Counter, defaultdict

import numpy as np

from cross_block_dual_r48_check import sample_order_ten
from envelope_block_cover_r27 import A6, A8, A9, projective_spins


def qnorm(a: np.ndarray) -> int:
    return max(abs(int(x @ a @ x)) for x in projective_spins(len(a)))


def audit(name: str, a: np.ndarray, m: int) -> None:
    n = len(a)
    q = qnorm(a)
    p = m / n
    p2 = m * (m - 1) / (n * (n - 1))
    p32q = p ** 1.5 * q
    words = []
    for x in projective_spins(n):
        e = int(x @ a @ x)
        for sigma in (-1, 1):
            words.append((sigma, x, sigma * e))

    by_qs = defaultdict(lambda: [0, 0])
    gaps = []
    automatic = 0
    favorable = 0
    selectors = 0
    for s in itertools.combinations(range(n), m):
        selectors += 1
        child = a[np.ix_(s, s)]
        qs = qnorm(child)
        ygap = qs - p32q
        gaps.append(ygap)
        if ygap <= 1e-12:
            automatic += 1
        local_fav = 0
        for sigma, x, oriented_full in words:
            local = sigma * int(x[list(s)] @ child @ x[list(s)])
            xcenter = local - p2 * oriented_full
            if ygap - xcenter <= 1e-12:
                local_fav += 1
        favorable += local_fav
        by_qs[qs][0] += local_fav
        by_qs[qs][1] += len(words)

    z0 = favorable / (selectors * len(words))
    print(
        name,
        "n", n,
        "m", m,
        "q", q,
        "Z0", round(z0, 9),
        "auto_selector_fraction", round(automatic / selectors, 9),
        "gap_range", (round(min(gaps), 6), round(max(gaps), 6)),
        "Q_S_hist", dict(sorted(Counter(round(g + p32q) for g in gaps).items())),
        "mass_by_Q_S", {
            qs: round(num / den, 9) for qs, (num, den) in sorted(by_qs.items())
        },
    )


def main() -> None:
    cases = [
        ("A6", A6, range(3, 6)),
        ("A8", A8, range(4, 8)),
        ("A9", A9, range(5, 9)),
        ("A10s0", sample_order_ten(0), range(5, 10)),
    ]
    for name, a, ms in cases:
        for m in ms:
            audit(name, a, m)


if __name__ == "__main__":
    main()
