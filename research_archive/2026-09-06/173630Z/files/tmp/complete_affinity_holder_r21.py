#!/usr/bin/env python3
"""Exact finite audits of complete-signing affinity aggregate identities."""

from __future__ import annotations

import itertools
import math

import numpy as np

from check_response_dual_r16 import A8, A9


def spins(n: int) -> np.ndarray:
    return np.array(list(itertools.product((-1, 1), repeat=n)), dtype=np.int64)


def lifted(A: np.ndarray):
    X = spins(len(A))
    base = np.einsum("bi,ij,bj->b", X, A, X)
    sigma = np.repeat(np.array((-1, 1), dtype=np.int64), len(X))
    XX = np.tile(X, (2, 1))
    energy = sigma * np.tile(base, 2)
    q = int(np.max(energy))
    return sigma, XX, energy, q


def deficit_sum(A: np.ndarray, beta: float) -> float:
    _, _, energy, q = lifted(A)
    return float(np.exp(-beta * (q - energy)).sum())


def audit(A: np.ndarray, name: str, beta: float) -> None:
    r = len(A)
    sigma, X, energy, q = lifted(A)
    deficit = q - energy
    weights = np.exp(-beta * deficit)
    D = float(weights.sum())
    nu = weights / D
    index = {(int(s), tuple(map(int, x))): j for j, (s, x) in enumerate(zip(sigma, X))}

    # Every complete edge has the exact deficit-at-most-four unfavorable
    # witness required by one-edge minimality.
    low = deficit <= 4
    for i in range(r):
        for j in range(i + 1, r):
            edge_sign = sigma * A[i, j] * X[:, i] * X[:, j]
            assert np.any(low & (edge_sign == -1))

    bcs = []
    kappas = []
    means = []
    child_product_rhs = []
    for i in range(r):
        h = sigma * X[:, i] * (X @ A[i])
        partners = []
        for s, x in zip(sigma, X):
            y = x.copy()
            y[i] *= -1
            partners.append(index[(int(s), tuple(map(int, y)))])
        partners = np.array(partners)
        assert np.array_equal(deficit[partners], deficit + 4 * h)
        bc_direct = float(np.exp(-beta * (deficit + 2 * h)).sum() / D)
        bc_sqrt = float(np.sqrt(weights * weights[partners]).sum() / D)
        assert np.isclose(bc_direct, bc_sqrt)
        bcs.append(bc_direct)
        kappas.append(-math.log(bc_direct) / beta)
        means.append(float(nu @ h))

        keep = [j for j in range(r) if j != i]
        C = A[np.ix_(keep, keep)]
        _, _, child_energy, qc = lifted(C)
        Dc = float(np.exp(-beta * (qc - child_energy)).sum())
        di = q - qc
        child_product_rhs.append(2 * math.exp(-beta * di) * Dc / D)
        assert np.isclose(bc_direct, child_product_rhs[-1])

    scaled_beta = beta * (1 - 2 / r)
    Dscaled = deficit_sum(A, scaled_beta)
    holder = math.exp(-2 * beta * q / r) * Dscaled / D
    geometric = math.exp(sum(math.log(v) for v in bcs) / r)
    arithmetic = sum(bcs) / r
    assert geometric + 1e-12 >= holder and arithmetic + 1e-12 >= holder
    assert np.isclose(np.prod(bcs), np.prod(child_product_rhs))
    assert np.isclose(sum(means), float(nu @ energy))

    expected_energy = float(nu @ energy)
    z = math.tanh(beta * (r - 1)) * expected_energy / (r * (r - 1))
    weak_lower = -math.log(max(1 - z, 1e-300)) / beta
    assert max(kappas) + 1e-10 >= weak_lower

    print(
        name,
        f"beta={beta}",
        {
            "Q": q,
            "BC_range": (min(bcs), max(bcs)),
            "kappa_range": (min(kappas), max(kappas)),
            "Holder": holder,
            "geo_BC": geometric,
            "mean_BC": arithmetic,
            "weak_max_kappa_lb": weak_lower,
        },
    )


def main() -> None:
    for beta in (0.1, 0.5, 1.0):
        audit(A8, "A8", beta)
        audit(A9, "A9", beta)
    print("PASS: complete energy, Holder, child-product, and pair-field audits")


if __name__ == "__main__":
    main()
