#!/usr/bin/env python3
"""Finite checks for the PSD compatibility-fibre gluing laws."""

from __future__ import annotations

from itertools import product
from math import sqrt

import numpy as np


def random_factor(rows: int, rank: int, rng: np.random.Generator) -> np.ndarray:
    Y = rng.normal(size=(rows, rank))
    Y *= sqrt(rows / max(float(np.sum(Y * Y)), 1e-30))
    return Y


def random_contraction(r1: int, r2: int, rng: np.random.Generator) -> np.ndarray:
    W = rng.normal(size=(r1, r2))
    norm = np.linalg.norm(W, 2)
    if norm > 1:
        W /= norm
    W *= rng.uniform(0.1, 1.0)
    return W


def q_metric(K: np.ndarray, L: np.ndarray) -> float:
    p = len(K)
    best = 0.0
    for eps0 in product((-1.0, 1.0), repeat=p):
        eps = np.asarray(eps0)
        best = max(best, abs(float(eps @ (K - L) @ eps)) / p**2)
    return best


def skinny_factor(K: np.ndarray, tolerance: float = 1e-9) -> np.ndarray:
    values, vectors = np.linalg.eigh(K)
    keep = values > tolerance
    return vectors[:, keep] * np.sqrt(values[keep])


def verify_douglas() -> int:
    rng = np.random.default_rng(281726)
    checks = 0
    for p1, p2, r1, r2 in ((3, 4, 2, 3), (5, 3, 1, 2), (4, 4, 3, 1)):
        for _ in range(20):
            Y1 = random_factor(p1, r1, rng)
            Y2 = random_factor(p2, r2, rng)
            W = random_contraction(r1, r2, rng)
            K1, K2 = Y1 @ Y1.T, Y2 @ Y2.T
            C = Y1 @ W @ Y2.T
            joined = np.block([[K1, C], [C.T, K2]])
            assert np.linalg.eigvalsh(joined).min() >= -1e-9

            recovered = np.linalg.pinv(Y1) @ C @ np.linalg.pinv(Y2).T
            assert np.allclose(recovered, W, atol=1e-8)
            assert np.linalg.norm(recovered, 2) <= 1 + 1e-8

            O1, _ = np.linalg.qr(rng.normal(size=(r1, r1)))
            O2, _ = np.linalg.qr(rng.normal(size=(r2, r2)))
            gauged = (Y1 @ O1) @ (O1.T @ W @ O2) @ (Y2 @ O2).T
            assert np.allclose(gauged, C, atol=1e-8)
            checks += 1
    return checks


def spectral_truncate(K: np.ndarray, threshold: float):
    values, vectors = np.linalg.eigh(K)
    keep = values > threshold
    P = vectors[:, keep] @ vectors[:, keep].T
    return P, P @ K @ P, int(keep.sum())


def verify_truncation() -> int:
    rng = np.random.default_rng(172628)
    checks = 0
    for p1, p2 in ((2, 3), (3, 3), (3, 4)):
        for tau in (0.05, 0.2, 0.6):
            for _ in range(20):
                r1, r2 = min(2, p1), min(2, p2)
                Y1 = random_factor(p1, r1, rng)
                Y2 = random_factor(p2, r2, rng)
                W = random_contraction(r1, r2, rng)
                K1, K2 = Y1 @ Y1.T, Y2 @ Y2.T
                C = Y1 @ W @ Y2.T
                K = np.block([[K1, C], [C.T, K2]])
                P1, H1, rank1 = spectral_truncate(K1, tau * p1)
                P2, H2, rank2 = spectral_truncate(K2, tau * p2)
                Ch = P1 @ C @ P2
                Kh = np.block([[H1, Ch], [Ch.T, H2]])
                assert np.linalg.eigvalsh(Kh).min() >= -1e-8
                assert rank1 <= int(np.floor(1 / tau))
                assert rank2 <= int(np.floor(1 / tau))
                exact = q_metric(K, Kh)
                assert exact <= sqrt(tau) + tau / 2 + 1e-8
                checks += 1
    return checks


def verify_square_root_example() -> int:
    checks = 0
    for s in (2, 3, 5):
        u = np.ones(s) / sqrt(s)
        for tau in (0.01, 0.09, 0.36):
            K1 = tau * s * np.outer(u, u)
            K2 = s * np.outer(u, u)
            C = sqrt(tau) * s * np.outer(u, u)
            K = np.block([[K1, C], [C.T, K2]])
            Kh = np.block(
                [
                    [np.zeros_like(K1), np.zeros_like(C)],
                    [np.zeros_like(C.T), K2],
                ]
            )
            assert np.linalg.eigvalsh(K).min() >= -1e-8
            expected = (tau + 2 * sqrt(tau)) / 4
            assert abs(q_metric(K, Kh) - expected) < 1e-10
            checks += 1
    return checks


def verify_factor_stability() -> int:
    rng = np.random.default_rng(71628)
    checks = 0
    for p1, p2, r1, r2 in ((2, 3, 2, 2), (3, 4, 2, 3)):
        for _ in range(100):
            Y1 = random_factor(p1, r1, rng)
            Y2 = random_factor(p2, r2, rng)
            W = random_contraction(r1, r2, rng)
            delta, zeta = rng.uniform(0.001, 0.08, size=2)

            E1 = rng.normal(size=Y1.shape)
            E1 *= delta * sqrt(p1) / max(np.linalg.norm(E1), 1e-30)
            E2 = rng.normal(size=Y2.shape)
            E2 *= delta * sqrt(p2) / max(np.linalg.norm(E2), 1e-30)
            Z1, Z2 = Y1 + E1, Y2 + E2
            # Rescale inward so the hypotheses ||Zi||_F<=sqrt(pi) hold.
            Z1 *= min(1.0, sqrt(p1) / np.linalg.norm(Z1))
            Z2 *= min(1.0, sqrt(p2) / np.linalg.norm(Z2))
            actual_delta = max(
                np.linalg.norm(Y1 - Z1) / sqrt(p1),
                np.linalg.norm(Y2 - Z2) / sqrt(p2),
            )

            E = rng.normal(size=W.shape)
            E *= zeta / max(np.linalg.norm(E), 1e-30)
            V = W + E
            norm = np.linalg.norm(V, 2)
            if norm > 1:
                V /= norm
            actual_zeta = np.linalg.norm(W - V, 2)

            K = np.block(
                [
                    [Y1 @ Y1.T, Y1 @ W @ Y2.T],
                    [Y2 @ W.T @ Y1.T, Y2 @ Y2.T],
                ]
            )
            L = np.block(
                [
                    [Z1 @ Z1.T, Z1 @ V @ Z2.T],
                    [Z2 @ V.T @ Z1.T, Z2 @ Z2.T],
                ]
            )
            assert np.linalg.eigvalsh(L).min() >= -1e-8
            assert q_metric(K, L) <= 2 * actual_delta + actual_zeta / 2 + 1e-8
            checks += 1
    return checks


def verify_multi_piece() -> int:
    rng = np.random.default_rng(62817)
    checks = 0
    ranks = (2, 1, 2)
    rows = (3, 2, 4)
    latent = 4
    # Gram blocks Ui Uj^T give a global PSD Omega with identity diagonals.
    isometries = []
    for rank in ranks:
        Q, _ = np.linalg.qr(rng.normal(size=(latent, rank)))
        isometries.append(Q.T)
    Omega = np.block(
        [[isometries[i] @ isometries[j].T for j in range(3)] for i in range(3)]
    )
    assert np.linalg.eigvalsh(Omega).min() >= -1e-9
    offset = 0
    for rank in ranks:
        block = Omega[offset : offset + rank, offset : offset + rank]
        assert np.allclose(block, np.eye(rank))
        offset += rank
    factors = [random_factor(row, rank, rng) for row, rank in zip(rows, ranks)]
    joined = np.block(
        [
            [factors[i] @ (isometries[i] @ isometries[j].T) @ factors[j].T for j in range(3)]
            for i in range(3)
        ]
    )
    assert np.linalg.eigvalsh(joined).min() >= -1e-8
    checks += 1
    return checks


def main() -> None:
    checks = 0
    checks += verify_douglas()
    checks += verify_truncation()
    checks += verify_square_root_example()
    checks += verify_factor_stability()
    checks += verify_multi_piece()
    print(f"PSD compatibility-fibre checks passed: {checks}")


if __name__ == "__main__":
    main()
