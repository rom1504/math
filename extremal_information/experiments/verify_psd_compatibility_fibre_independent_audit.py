#!/usr/bin/env python3
"""Independent checks for the PSD compatibility-fibre draft.

Adds converse recovery of a multi-piece correlation carrier, a zero-rank
Douglas edge case, and the tau=1 endpoint diagnostic for PF.21--PF.22.
"""

from __future__ import annotations

from math import ceil, log, sqrt

import numpy as np


def skinny_factor(matrix: np.ndarray, tolerance: float = 1e-10) -> np.ndarray:
    values, vectors = np.linalg.eigh(matrix)
    keep = values > tolerance
    return vectors[:, keep] * np.sqrt(values[keep])


def q_metric(left: np.ndarray, right: np.ndarray) -> float:
    p = len(left)
    best = 0.0
    for mask in range(1 << p):
        word = np.asarray([1 if (mask >> i) & 1 else -1 for i in range(p)])
        best = max(best, abs(float(word @ (left - right) @ word)) / p**2)
    return best


def check_arbitrary_multi_piece_recovery() -> int:
    rng = np.random.default_rng(170826)
    checks = 0
    for ranks, rows, latent in (
        ((1, 2, 1), (3, 4, 2), 4),
        ((2, 2, 3), (4, 3, 5), 6),
    ):
        local_factors = []
        row_isometries = []
        global_blocks = []
        for rank, row in zip(ranks, rows):
            y = rng.normal(size=(row, rank))
            while np.linalg.matrix_rank(y) < rank:
                y = rng.normal(size=(row, rank))
            q, _ = np.linalg.qr(rng.normal(size=(latent, rank)))
            u = q.T  # u u^T = I_rank.
            local_factors.append(y)
            row_isometries.append(u)
            global_blocks.append(y @ u)

        global_factor = np.vstack(global_blocks)
        joined = global_factor @ global_factor.T
        recovered_local = []
        offset = 0
        for row in rows:
            block = joined[offset : offset + row, offset : offset + row]
            recovered_local.append(skinny_factor(block))
            offset += row

        block_diagonal = np.zeros((sum(rows), sum(ranks)))
        row_offset = 0
        rank_offset = 0
        for row, rank, factor in zip(rows, ranks, recovered_local):
            block_diagonal[
                row_offset : row_offset + row,
                rank_offset : rank_offset + rank,
            ] = factor
            row_offset += row
            rank_offset += rank

        inverse = np.linalg.pinv(block_diagonal)
        omega = inverse @ joined @ inverse.T
        assert np.linalg.eigvalsh(omega).min() >= -1e-8
        rank_offset = 0
        for rank in ranks:
            diagonal = omega[
                rank_offset : rank_offset + rank,
                rank_offset : rank_offset + rank,
            ]
            assert np.allclose(diagonal, np.eye(rank), atol=1e-8)
            rank_offset += rank
        assert np.allclose(block_diagonal @ omega @ block_diagonal.T, joined, atol=1e-8)
        checks += 1
    return checks


def check_zero_rank_douglas() -> int:
    rng = np.random.default_rng(281726)
    k2_factor = rng.normal(size=(4, 2))
    k2 = k2_factor @ k2_factor.T
    k1 = np.zeros((3, 3))
    cross = np.zeros((3, 4))
    joined = np.block([[k1, cross], [cross.T, k2]])
    assert np.linalg.eigvalsh(joined).min() >= -1e-10
    y1 = skinny_factor(k1)
    y2 = skinny_factor(k2)
    assert y1.shape == (3, 0)
    contraction = np.zeros((0, y2.shape[1]))
    assert np.array_equal(y1 @ contraction @ y2.T, cross)

    # Any nonzero cross block violates PSD when the first marginal is zero.
    bad_cross = cross.copy()
    bad_cross[0, 0] = 1e-3
    bad = np.block([[k1, bad_cross], [bad_cross.T, k2]])
    assert np.linalg.eigvalsh(bad).min() < -1e-10
    return 2


def check_square_root_endpoint() -> int:
    s = 3
    u = np.ones(s) / sqrt(s)
    tau = 1.0
    k1 = tau * s * np.outer(u, u)
    k2 = s * np.outer(u, u)
    cross = sqrt(tau) * s * np.outer(u, u)
    joined = np.block([[k1, cross], [cross.T, k2]])

    # At strict cutoff tau*p_i both equal-threshold eigenvalues are removed.
    truncated = np.zeros_like(joined)
    observed = q_metric(joined, truncated)
    claimed_pf22_expression = (tau + 2 * sqrt(tau)) / 4
    assert abs(observed - 1.0) < 1e-12
    assert abs(claimed_pf22_expression - 0.75) < 1e-12
    assert observed != claimed_pf22_expression
    return 1


def check_net_and_response_arithmetic() -> int:
    checks = 0
    for eta in (0.05, 0.2, 1.0):
        tau = (eta / 4) ** 2
        delta = zeta = eta / 16
        metric_error = 2 * sqrt(tau) + tau + 4 * delta + zeta
        assert metric_error <= 7 * eta / 8 + 1e-15
        rank = ceil(1 / tau)
        assert rank >= 16 / eta**2
        radius_ratio = 2 * sqrt(16 / eta**2) / zeta
        assert abs(radius_ratio - 128 / eta**2) < 1e-10
        log_bound = 512 / eta**4 * log(1 + 128 / eta**2)
        assert log_bound > 0

        c = 1.0
        response_error = c * sqrt(metric_error / 2) + c * c * metric_error / 8
        assert response_error > 0
        checks += 1
    return checks


def main() -> None:
    multi = check_arbitrary_multi_piece_recovery()
    zero = check_zero_rank_douglas()
    endpoint = check_square_root_endpoint()
    arithmetic = check_net_and_response_arithmetic()
    print(
        "PSD compatibility-fibre independent audit: PASS",
        f"multi_recovery={multi}",
        f"zero_rank={zero}",
        f"tau_endpoint_diagnostic={endpoint}",
        f"arithmetic={arithmetic}",
    )


if __name__ == "__main__":
    main()
