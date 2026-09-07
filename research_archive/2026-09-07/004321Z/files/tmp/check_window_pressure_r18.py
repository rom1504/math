#!/usr/bin/env python3
"""Exact-energy/numerical-pressure checks for the Wave-18 window telescope."""

from __future__ import annotations

import math

import numpy as np

from check_selected_child_soft_r17 import A9, oriented_energies, p0


def child(a: np.ndarray, i: int) -> np.ndarray:
    ids = [j for j in range(len(a)) if j != i]
    return a[np.ix_(ids, ids)]


def qnorm(a: np.ndarray) -> int:
    e = oriented_energies(a)
    return int(e.max())


def deficit_sum(a: np.ndarray, beta: float) -> float:
    e = oriented_energies(a).astype(float)
    q = float(e.max())
    return float(np.exp(-beta * (q - e)).sum())


def audit(beta: float) -> None:
    # The named flat chain A9 -> A9[-7] -> A9[-{7,0}].
    mats = [A9, child(A9, 7)]
    mats.append(child(mats[-1], 0))
    alpha = 24 / (9 ** 1.5)
    root_p = p0(mats[0], beta)
    print(f"beta={beta:.12g}")
    for a in mats:
        m = len(a)
        q = qnorm(a)
        p = p0(a, beta)
        reward = root_p - p
        target = alpha * (9 ** 1.5 - m ** 1.5)
        shortfall = target - reward
        entropy = (m * math.log(2) - math.log(deficit_sum(a, beta))) / beta
        assert abs(q - (p + entropy)) < 2e-10
        exact_excess = q - alpha * m ** 1.5
        assert abs(exact_excess - ((p - alpha * m ** 1.5) + entropy)) < 2e-10
        print(
            f"  m={m} Q={q:2d} P={p: .9f} reward={reward: .9f} "
            f"pressure-shortfall={shortfall: .9f} entropy={entropy: .9f} "
            f"exact-excess={exact_excess: .9f}"
        )


def main() -> None:
    for b in (0.5, 1, 2, 4, 8, 16):
        audit(b / 3)  # beta=b/sqrt(9), frozen over the whole window.
    print("PASS: fixed-temperature telescopes and entropy corrections agree")


if __name__ == "__main__":
    main()
