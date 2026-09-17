#!/usr/bin/env python3
"""Exact finite diagnostics for Wave 37 edge-flip selector compression."""

from __future__ import annotations

import itertools
import sys

import numpy as np
from scipy.optimize import linprog

sys.path.insert(0, "/home/math/quadra/tmp")
from envelope_block_cover_r27 import A6, A8, A9


def projective_spins(n: int):
    for tail in itertools.product((-1, 1), repeat=n - 1):
        yield np.asarray((1,) + tail, dtype=np.int64)


def q_value(a: np.ndarray) -> int:
    return max(abs(int(x @ a @ x)) for x in projective_spins(len(a)))


def fractional_cover(inc: np.ndarray) -> float | None:
    if np.any(inc.sum(axis=1) == 0):
        return None
    sol = linprog(
        np.ones(inc.shape[1]),
        A_ub=-inc,
        b_ub=-np.ones(inc.shape[0]),
        bounds=(0.0, None),
        method="highs",
    )
    assert sol.success
    return float(sol.fun)


def audit(a: np.ndarray, m: int, name: str, expected_tau: float,
          expected_mean_cap_misses: int) -> None:
    n = len(a)
    q = q_value(a)
    p = m / n
    p2 = m * (m - 1) / (n * (n - 1))
    allowance = (p ** 1.5 - p2) * q
    assert p2 >= 0.5

    states = []
    for x in projective_spins(n):
        row = int((a @ x) @ (a @ x))
        for sigma in (-1, 1):
            energy = sigma * int(x @ a @ x)
            states.append((x, sigma, energy, q - energy, row))

    selectors = list(itertools.combinations(range(n), m))
    incidence = np.zeros((len(selectors), len(states)), dtype=np.int8)
    favorable = np.zeros_like(incidence)

    for si, s in enumerate(selectors):
        idx = np.asarray(s, dtype=int)
        child = a[np.ix_(idx, idx)]
        q_child = q_value(child)
        assert q_child <= q
        for di, (x, sigma, energy, deficit, _) in enumerate(states):
            y = x[idx]
            child_energy = sigma * int(y @ child @ y)
            # Flipping every edge outside E(S) changes E to 2*c_S-E.
            perturbed_energy = 2 * child_energy - energy
            certifies = perturbed_energy >= q
            incidence[si, di] = certifies
            h = q_child - child_energy - p2 * deficit - allowance
            favorable[si, di] = h <= 1e-10
            if certifies:
                assert child_energy + 1e-12 >= q - deficit / 2
                assert h <= -allowance + 1e-10

    assert np.all(incidence.sum(axis=1) >= 1)
    assert np.all(incidence <= favorable)

    all_tau = fractional_cover(incidence)
    assert all_tau is not None and abs(all_tau - expected_tau) < 1e-8

    mean_cap = n * (n - 1)
    good = np.asarray([row <= mean_cap for *_, row in states])
    misses = int(np.sum(incidence[:, good].sum(axis=1) == 0))
    assert misses == expected_mean_cap_misses

    print(
        name,
        "n", n,
        "m", m,
        "p2", p2,
        "selectors", len(selectors),
        "states", len(states),
        "witness-degree", (int(incidence.sum(1).min()), int(incidence.sum(1).max())),
        "fractional-cover", all_tau,
        "mean-row-cap-misses", misses,
    )


def main() -> None:
    audit(A6, 5, "A6", 2.0, 0)
    audit(A8, 6, "A8", 10.0, 4)
    audit(A9, 7, "A9", 13.75, 5)
    print("PASS discrete_minimality_compression_r37_check")


if __name__ == "__main__":
    main()
