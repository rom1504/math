#!/usr/bin/env python3
"""Checks the finite identities in common_capacity_r34.md.

The signing geometry is exhaustively enumerated.  Floating point is used only
for exponentials and p^(3/2); all energy, distance, and row-square data are
integer.
"""

from __future__ import annotations

import itertools
import math
import sys

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from envelope_block_cover_r27 import A6, A8, A9, projective_spins


def spins(n: int) -> np.ndarray:
    return np.array(list(projective_spins(n)), dtype=np.int64)


def projective_distance(x: np.ndarray, y: np.ndarray) -> int:
    h = int(np.count_nonzero(x != y))
    return min(h, len(x) - h)


def audit(A: np.ndarray, m: int) -> None:
    n = len(A)
    xs = spins(n)
    full_energy = np.einsum("bi,ij,bj->b", xs, A, xs)
    qn = int(np.max(np.abs(full_energy)))
    row = np.einsum("bi,ij,bj->b", xs, A @ A, xs)
    centers = xs[row <= 2 * n * (n - 1)]
    assert len(centers) > 0

    selectors = list(itertools.combinations(range(n), m))
    p = m / n
    p2 = m * (m - 1) / (n * (n - 1))
    B = (p ** 1.5 - p2) * qn
    tolerance = 0.0
    ys = spins(m)

    qS = []
    distances = np.empty((len(centers), len(selectors)), dtype=np.int64)
    delta_abs = np.empty_like(distances)
    for j, selector in enumerate(selectors):
        S = list(selector)
        AS = A[np.ix_(S, S)]
        child_energy = np.einsum("bi,ij,bj->b", ys, AS, ys)
        q_child = int(np.max(np.abs(child_energy)))
        qS.append(q_child)
        favorable = ys[(q_child - np.abs(child_energy)) <= B + tolerance + 1e-12]
        assert len(favorable)
        for i, z in enumerate(centers):
            ez = int(z[S] @ AS @ z[S])
            delta_abs[i, j] = q_child - abs(ez)
            distances[i, j] = min(projective_distance(z[S], y) for y in favorable)

    # B34.13, including the exact zero-distance implication B34.15.
    assert np.all(
        np.maximum(delta_abs - B - tolerance, 0.0)
        <= 4 * (m - 1) * distances + 1e-12
    )
    assert np.all(distances[delta_abs <= B + tolerance + 1e-12] == 0)

    lam = 0.73
    s = 3
    K = np.exp(-lam * distances)
    qz = K.mean(axis=1)
    L = np.mean(qz ** s)
    qstar = float(np.max(qz))
    N = len(centers)
    assert qstar ** s / N <= L + 1e-14
    assert L <= qstar ** s + 1e-14

    # Renyi identity B34.8.
    Z = float(np.mean(qz))
    density = qz / Z
    Ds = math.log(float(np.mean(density ** s))) / (s - 1)
    rhs = Z ** s * math.exp((s - 1) * Ds)
    assert abs(L - rhs) <= 2e-13

    # Per-center hard/soft sandwich B34.9--B34.10.
    for d in range(int(np.max(distances)) + 1):
        gd = np.mean(distances <= d, axis=1)
        lower = np.exp(-lam * d) * gd
        upper = gd + np.exp(-lam * (d + 1)) * (1 - gd)
        assert np.all(lower <= qz + 2e-14)
        assert np.all(qz <= upper + 2e-14)

    # Positive-part deficit obstruction B34.14.
    gamma = lam / (4 * (m - 1))
    deficit_laplace = np.mean(
        np.exp(-gamma * np.maximum(delta_abs - B - tolerance, 0.0)), axis=1
    )
    assert np.all(qz <= deficit_laplace + 2e-14)

    # Parent-slack identities B34.16, for both orientations and every state.
    qbar = float(np.mean(qS))
    G = qbar - p ** 1.5 * qn
    for z, energy in zip(xs, full_energy):
        child_raw = []
        for selector in selectors:
            S = list(selector)
            child_raw.append(int(z[S] @ A[np.ix_(S, S)] @ z[S]))
        child_raw = np.array(child_raw)
        for sigma in (-1, 1):
            Delta = qn - sigma * int(energy)
            delta = np.array(qS, dtype=float) - sigma * child_raw
            assert abs(float(np.mean(delta)) - (B + G + p2 * Delta)) < 2e-12
            effective = delta - p2 * Delta - B
            assert abs(float(np.mean(effective)) - G) < 2e-12

    print(
        f"PASS n={n} m={m} centers={N} selectors={len(selectors)} "
        f"q*={qstar:.9f} L={L:.9f} Z={Z:.9f} D_s={Ds:.9f}"
    )


def synthetic_rare_center() -> None:
    # Equality case showing the Renyi recovery of all but one rarity factor.
    N, s, value = 37, 5, 0.4
    q = np.zeros(N)
    q[0] = value
    Z = float(np.mean(q))
    L = float(np.mean(q ** s))
    density = q / Z
    Ds = math.log(float(np.mean(density ** s))) / (s - 1)
    assert abs(Ds - math.log(N)) < 1e-13
    assert abs(L - value ** s / N) < 1e-15
    assert abs(L - Z ** s * math.exp((s - 1) * Ds)) < 1e-14
    print("PASS synthetic one-center Renyi extremizer")


if __name__ == "__main__":
    synthetic_rare_center()
    audit(A6, 3)
    audit(A8, 4)
    audit(A9, 5)
    print("ALL COMMON-CAPACITY CHECKS PASS")
