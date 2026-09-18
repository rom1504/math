#!/usr/bin/env python3
"""Wave 24 audits: full competitor Fourier law and plaquette affinity bound."""

from __future__ import annotations

import itertools
import math

import numpy as np

from check_finite_bridge_r16 import A5
from check_response_dual_r16 import A8, A9, spins


def oriented_projective(A: np.ndarray):
    n = len(A)
    X0 = spins(n)
    X = np.tile(X0, (2, 1))
    sigma = np.repeat(np.array((-1, 1), dtype=np.int64), len(X0))
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    signs = np.column_stack(
        [sigma * A[i, j] * X[:, i] * X[:, j] for i, j in edges]
    )
    energy = 2 * signs.sum(axis=1)
    q = int(energy.max())
    deficit = q - energy
    return X, sigma, edges, signs, deficit, q


def fwht(values: np.ndarray) -> np.ndarray:
    out = np.asarray(values, dtype=float).copy()
    step = 1
    while step < len(out):
        for start in range(0, len(out), 2 * step):
            left = out[start : start + step].copy()
            right = out[start + step : start + 2 * step].copy()
            out[start : start + step] = left + right
            out[start + step : start + 2 * step] = left - right
        step *= 2
    return out


def fourier_audit(beta: float) -> None:
    A = A5
    n = len(A)
    _, _, edges, signs, deficit, q = oriented_projective(A)
    N = len(edges)
    Dvals = np.empty(1 << N)
    for mask in range(1 << N):
        chosen = [k for k in range(N) if (mask >> k) & 1]
        shift = 4 * signs[:, chosen].sum(axis=1) if chosen else 0
        Dvals[mask] = np.exp(-beta * (deficit + shift)).sum()
    # A5 is an exact minimizer, hence every perturbation sum has a term >=1.
    assert Dvals.min() >= 1 - 1e-11

    transform = fwht(Dvals) / (1 << N)
    c = 2**n * math.exp(-beta * q)
    predicted = np.zeros_like(transform)
    for mask in range(1 << N):
        k = bin(mask).count("1")
        degree = [0] * n
        edge_product = 1
        for eidx, (i, j) in enumerate(edges):
            if (mask >> eidx) & 1:
                degree[i] ^= 1
                degree[j] ^= 1
                edge_product *= int(A[i, j])
        if k % 2 == 0 and not any(degree):
            predicted[mask] = (
                c
                * math.cosh(2 * beta) ** (N - k)
                * math.sinh(2 * beta) ** k
                * edge_product
            )
    assert np.allclose(transform, predicted, rtol=2e-9, atol=2e-9)
    assert math.isclose(
        Dvals.mean(),
        2**n * math.exp(-beta * q) * math.cosh(2 * beta) ** N,
        rel_tol=2e-11,
    )
    assert math.isclose(
        np.mean(Dvals**2), np.sum(transform**2), rel_tol=2e-11
    )
    print(
        "A5 full-edge Fourier",
        f"beta={beta}",
        {
            "N": N,
            "Q": q,
            "min_D_T": float(Dvals.min()),
            "max_D_T": float(Dvals.max()),
            "nonzero_Fourier_coefficients": int(np.count_nonzero(np.abs(transform) > 1e-9)),
            "even_Eulerian_support": int(np.count_nonzero(predicted)),
        },
    )


def affinity_data(A: np.ndarray, beta: float):
    n = len(A)
    _, _, edges, signs, deficit, q = oriented_projective(A)
    weights = np.exp(-beta * deficit)
    D = float(weights.sum())
    bcs = []
    for i in range(n):
        H = [k for k, (u, v) in enumerate(edges) if i in (u, v)]
        h = signs[:, H].sum(axis=1)
        bcs.append(float(np.exp(-beta * (deficit + 2 * h)).sum() / D))
    return q, D, np.asarray(bcs)


def plaquette_constant(beta: float) -> float:
    lam = math.exp(4 * beta)
    if lam <= 3:
        return 4 * math.sqrt(lam) / (lam + 1)
    return math.sqrt(2 * lam / (lam - 1))


def plaquette_audit(A: np.ndarray, name: str, beta: float) -> None:
    n = len(A)
    q, D_projective, bcs = affinity_data(A, beta)
    Cbeta = plaquette_constant(beta)
    assert max(bcs[i] + bcs[j] for i in range(n) for j in range(i + 1, n)) <= Cbeta + 1e-11

    # Exhaust every conditional two-spin square, including both orientations.
    max_local = 0.0
    for sigma in (-1, 1):
        for i in range(n):
            for j in range(i + 1, n):
                rest = [k for k in range(n) if k not in (i, j)]
                for tail in itertools.product((-1, 1), repeat=n - 2):
                    f = {}
                    for xi, xj in itertools.product((-1, 1), repeat=2):
                        x = np.ones(n, dtype=np.int64)
                        x[rest] = tail
                        x[i], x[j] = xi, xj
                        energy = sigma * int(x @ A @ x)
                        f[(xi, xj)] = math.exp(beta * energy / 2)
                    a = f[(-1, -1)]
                    b = f[(1, -1)]
                    c = f[(-1, 1)]
                    d = f[(1, 1)]
                    cross = a * d / (b * c)
                    assert math.isclose(
                        abs(math.log(cross)), 4 * beta, rel_tol=2e-11, abs_tol=2e-11
                    )
                    z = a * a + b * b + c * c + d * d
                    local = 2 * ((a * b + c * d) + (a * c + b * d)) / z
                    max_local = max(max_local, local)
                    assert local <= Cbeta + 2e-11

    # Signed even-cycle polynomial and its unsigned majorant, computed from
    # the primal partition sums rather than enumerating a huge dual code.
    all_plus = np.ones((n, n), dtype=np.int64) - np.eye(n, dtype=np.int64)
    # oriented_projective used its own q; reconstruct weights with baseline q.
    X0 = spins(n)
    raw_plus = np.einsum("bi,ij,bj->b", X0, all_plus, X0)
    plus_energy = np.r_[ -raw_plus, raw_plus ]
    D_plus_same_q = float(np.exp(-beta * (q - plus_energy)).sum())
    N = n * (n - 1) // 2
    base = 2**n * math.exp(-beta * q) * math.cosh(2 * beta) ** N
    signed_cycle = D_projective / base
    unsigned_cycle = D_plus_same_q / base
    assert abs(signed_cycle) <= unsigned_cycle + 1e-10

    print(
        name,
        f"beta={beta}",
        {
            "Q": q,
            "BC_range": (float(bcs.min()), float(bcs.max())),
            "max_pair_BC_sum": float(
                max(bcs[i] + bcs[j] for i in range(n) for j in range(i + 1, n))
            ),
            "plaquette_Cbeta": Cbeta,
            "universal_max_kappa_lb": math.log(2 / Cbeta) / beta,
            "max_conditional_square_value": max_local,
            "signed_cycle_polynomial": signed_cycle,
            "unsigned_cycle_majorant": unsigned_cycle,
            "cancellation_ratio": signed_cycle / unsigned_cycle,
        },
    )


def main() -> None:
    for beta in (0.5, 1.0):
        fourier_audit(beta)
        plaquette_audit(A8, "A8", beta)
        plaquette_audit(A9, "A9", beta)
    print("PASS: full competitor Fourier law and complete-signing plaquette bound")


if __name__ == "__main__":
    main()
