#!/usr/bin/env python3
"""Audit edge-witness incidence and coordinate affinities."""

from __future__ import annotations

import itertools
import math

import numpy as np
from scipy.optimize import linprog

from check_response_dual_r16 import A9, spins


def abstract_ball_affinity(r: int, k: int) -> float:
    numerator = 2 * sum(math.comb(r - 1, j) for j in range(k))
    denominator = sum(math.comb(r, j) for j in range(k + 1))
    return numerator / denominator


def audit_abstract() -> None:
    for r in (25, 64, 121, 256):
        k = int(math.isqrt(r))
        bc = abstract_ball_affinity(r, k)
        # The negative-orientation all-positive center has s_e=-1 on every
        # edge, so one deficit-zero abstract witness covers the full edge set.
        print(f"abstract r={r} k={k} BC={bc:.9f} reward/log-scale={-math.log(bc):.6f}")
        assert bc > math.exp(-math.sqrt(r))


def all_spins(n: int) -> np.ndarray:
    return np.array(list(itertools.product((-1, 1), repeat=n)), dtype=np.int64)


def audit_a9() -> None:
    A = A9
    n = len(A)
    # Exact projective incidence LP.
    X = spins(n)
    base_energy = np.einsum("bi,ij,bj->b", X, A, X)
    q = int(np.max(np.abs(base_energy)))
    sigmas = np.repeat(np.array((-1, 1), dtype=np.int64), len(X))
    XX = np.tile(X, (2, 1))
    energies = sigmas * np.tile(base_energy, 2)
    deficits = q - energies

    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    low = np.flatnonzero(deficits <= 4)
    incidence = np.zeros((len(edges), len(low)), dtype=float)
    for eidx, (i, j) in enumerate(edges):
        signs = sigmas[low] * A[i, j] * XX[low, i] * XX[low, j]
        incidence[eidx] = signs == -1
        assert incidence[eidx].any()

    # Exact primal certificate: listed weights are divided by three.
    primal_indices = (0, 4, 5, 6, 10, 15, 19, 20, 22)
    primal_numerators = np.array((1, 2, 1, 1, 2, 1, 2, 1, 1), dtype=int)
    assert primal_numerators.sum() == 12
    assert np.all(incidence[:, primal_indices] @ primal_numerators >= 3)
    # Exact dual certificate: every low state is negative on at most one of
    # these four edges, each assigned dual weight one.
    dual_edges = ((0, 2), (0, 4), (2, 5), (4, 5))
    dual_indices = [edges.index(edge) for edge in dual_edges]
    assert np.all(incidence[dual_indices].sum(axis=0) == 1)

    result = linprog(
        np.ones(len(low)),
        A_ub=-incidence,
        b_ub=-np.ones(len(edges)),
        bounds=(0, None),
        method="highs",
    )
    assert result.success
    print(
        "A9",
        {"Q": q, "projective_low_states": len(low), "fractional_cover_exact": 4,
         "min_edge_multiplicity": int(incidence.sum(axis=1).min()),
         "max_edge_multiplicity": int(incidence.sum(axis=1).max())},
    )

    # Full lifted states for the Bhattacharyya ratio.  This duplicates every
    # oriented-projective state twice and leaves (10.623) unchanged.
    X = all_spins(n)
    base_energy = np.einsum("bi,ij,bj->b", X, A, X)
    sigmas = np.repeat(np.array((-1, 1), dtype=np.int64), len(X))
    XX = np.tile(X, (2, 1))
    energies = sigmas * np.tile(base_energy, 2)
    deficits = q - energies
    index = {(int(sigma), tuple(map(int, x))): idx for idx, (sigma, x) in enumerate(zip(sigmas, XX))}
    for beta in (0.1, 0.5, 1.0):
        weights = np.exp(-beta * deficits)
        denom = float(weights.sum())
        kappas = []
        for i in range(n):
            partners = []
            for sigma, x in zip(sigmas, XX):
                y = x.copy()
                y[i] *= -1
                partners.append(index[(int(sigma), tuple(map(int, y)))])
            bc = float(np.sqrt(weights * weights[np.array(partners)]).sum() / denom)
            kappas.append(-math.log(bc) / beta)
        assert max(kappas) < 3
        print(f"A9 beta={beta}: kappa range=({min(kappas):.9f},{max(kappas):.9f})")


def main() -> None:
    audit_abstract()
    audit_a9()
    print("PASS: abstract affinity and actual A9 witness/affinity audits")


if __name__ == "__main__":
    main()
