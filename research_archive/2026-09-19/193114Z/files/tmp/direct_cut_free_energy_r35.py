#!/usr/bin/env python3
"""Exact/numerical checks for the Wave 35 direct-cut hinge free energy.

All cut, selector, energy, and row-square data are exhaustively enumerated.
Floating point enters only through p^(3/2), exponentials, and logarithms.
"""

from __future__ import annotations

import itertools
import math
import sys

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from envelope_block_cover_r27 import A6, A8, A9, all_data


def canonical(x: np.ndarray) -> tuple[int, ...]:
    x = x.copy()
    if x[0] < 0:
        x *= -1
    return tuple(int(v) for v in x)


def full_arrays(A: np.ndarray, m: int):
    n = len(A)
    q, selectors, records, q_child, den = all_data(A, m)
    rows = np.array([r[2] for r in records], dtype=np.int64)
    xnum = np.stack([r[3] for r in records], axis=0)
    y = q_child.astype(float) - q * (m / n) ** 1.5
    effective = y[None, :] - xnum / den
    keys = {(r[0], r[1]): i for i, r in enumerate(records)}
    return q, selectors, records, rows, effective, keys


def audit_block_increments(A: np.ndarray, m: int) -> None:
    n = len(A)
    q, selectors, records, rows, effective, keys = full_arrays(A, m)
    p2 = m * (m - 1) / (n * (n - 1))
    rng = np.random.default_rng(3504 + n)
    for _ in range(300):
        j = int(rng.integers(len(records)))
        x_tuple, sigma, row, _ = records[j]
        x = np.array(x_tuple, dtype=np.int64)
        H = set(np.flatnonzero(rng.integers(0, 2, n)).tolist())
        xp = x.copy()
        if H:
            xp[list(H)] *= -1
        jp = keys[(canonical(xp), sigma)]

        shore = sum(
            int(A[i, k] * x[i] * x[k])
            for i in H
            for k in range(n)
            if k not in H
        )
        v = A[:, list(H)] @ x[list(H)] if H else np.zeros(n, dtype=np.int64)
        delta_row = (
            -4 * sum(int(x[i] * (A @ A @ x)[i]) for i in H)
            + 4 * int(v @ v)
        )
        assert rows[jp] - row == delta_row

        for k, S0 in enumerate(selectors):
            S = set(S0)
            child_shore = sum(
                int(A[i, ell] * x[i] * x[ell])
                for i in H & S
                for ell in S
                if ell not in H
            )
            delta = 4 * sigma * (child_shore - p2 * shore)
            assert abs((effective[jp, k] - effective[j, k]) - delta) < 2e-12

        jo = keys[(x_tuple, -sigma)]
        # If X=Y-effective, orientation reversal changes effective by 2X.
        y = np.array([
            max(
                abs(int(z @ A[np.ix_(S, S)] @ z))
                for z in itertools.product((-1, 1), repeat=m)
                for z in [np.array(z, dtype=np.int64)]
            )
            - q * (m / n) ** 1.5
            for S in selectors
        ])
        assert np.max(np.abs(effective[jo] - effective[j] - 2 * (y - effective[j]))) < 2e-12


def audit_free_energy(A: np.ndarray, name: str, m: int) -> None:
    n = len(A)
    _, _, records, rows, effective, keys = full_arrays(A, m)
    threshold = 0.0
    phi = np.maximum(effective - threshold, 0.0)
    hard = np.mean(effective <= threshold + 1e-12, axis=1)

    for lam in (0.05, 0.2, 1.0, 5.0):
        kernel = np.exp(-lam * phi)
        z = kernel.mean(axis=1)
        for gamma in (0.0, 0.1, 1.0):
            objective = -np.log(z) + gamma * rows / (n * (n - 1))
            j = int(np.argmin(objective))
            x_tuple, sigma, _, _ = records[j]
            x = np.array(x_tuple, dtype=np.int64)
            posterior = kernel[j] / kernel[j].sum()

            # Gibbs variational equality for the selected cut.
            kl = float(np.sum(posterior * np.log(posterior * len(posterior))))
            assert abs(math.log(z[j]) + lam * float(posterior @ phi[j]) + kl) < 2e-12

            # Global optimality gives every physical-vertex and orientation
            # Euler inequality.  Re-gauging x after a physical flip does not
            # change its projective cut.
            neighbors = []
            for i in range(n):
                xp = x.copy()
                xp[i] *= -1
                neighbors.append(keys[(canonical(xp), sigma)])
            neighbors.append(keys[(x_tuple, -sigma)])
            for jp in neighbors:
                eta = phi[jp] - phi[j]
                ratio = float(posterior @ np.exp(-lam * eta))
                delta_row = int(rows[jp] - rows[j])
                assert math.log(ratio) <= gamma * delta_row / (n * (n - 1)) + 2e-12

            # Continuous-threshold soft/hard sandwich for several buffers.
            for u in (0.5, 2.0, 5.0):
                enlarged = float(np.mean(effective[j] <= threshold + u + 1e-12))
                assert z[j] <= enlarged + math.exp(-lam * u) + 2e-12

            print({
                "A": name,
                "m": m,
                "lambda": lam,
                "gamma": gamma,
                "row": int(rows[j]),
                "hard": round(float(hard[j]), 6),
                "hard_best": round(float(hard.max()), 6),
                "soft": round(float(z[j]), 6),
                "KL": round(kl, 6),
            })


if __name__ == "__main__":
    for A, name, m in ((A6, "A6", 3), (A8, "A8", 4), (A9, "A9", 5)):
        audit_block_increments(A, m)
        audit_free_energy(A, name, m)
    print("PASS direct_cut_free_energy_r35")
