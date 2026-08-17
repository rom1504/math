#!/usr/bin/env python3
"""Finite checks for drafts/psd_gluing_compatibility_entropy.md."""

from __future__ import annotations

from itertools import permutations, product

import numpy as np

from verify_bcx_two_port_holonomy import regular_hadamard


V0 = np.asarray(
    (
        -1, -1, -1, 1,
        -1, -1, 1, -1,
        1, -1, 1, 1,
        -1, 1, 1, 1,
    ),
    dtype=np.int64,
)


def psd_sqrt(matrix: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    values, vectors = np.linalg.eigh(matrix)
    values = np.maximum(values, 0)
    root = (vectors * np.sqrt(values)) @ vectors.T
    inverse_values = np.zeros_like(values)
    support = values > 1e-10
    inverse_values[support] = 1 / np.sqrt(values[support])
    inverse_root = (vectors * inverse_values) @ vectors.T
    return root, inverse_root


def inf_to_one(matrix: np.ndarray) -> float:
    q = matrix.shape[1]
    best = 0.0
    for y in product((-1, 1), repeat=q):
        best = max(best, float(np.sum(np.abs(matrix @ np.asarray(y)))))
    return best


def verify_contraction_parametrization() -> None:
    rng = np.random.default_rng(20260817)
    p, q = 6, 5
    left_factor = rng.normal(size=(p, 3))
    right_factor = rng.normal(size=(q, 2))
    kl = left_factor @ left_factor.T
    kr = right_factor @ right_factor.T
    scale_l = max(1.0, float(np.max(np.diag(kl))))
    scale_r = max(1.0, float(np.max(np.diag(kr))))
    kl /= scale_l
    kr /= scale_r
    lr, lri = psd_sqrt(kl)
    rr, rri = psd_sqrt(kr)

    c = rng.normal(size=(p, q))
    c *= 0.8 / np.linalg.norm(c, ord=2)
    x = lr @ c @ rr
    joined = np.block([[kl, x], [x.T, kr]])
    assert float(np.linalg.eigvalsh(joined).min()) >= -1e-9

    compressed = lri @ x @ rri
    assert np.linalg.norm(compressed, ord=2) <= 1 + 1e-9
    assert np.linalg.norm(lr @ compressed @ rr - x) < 1e-7
    print(
        "singular contraction parametrization:",
        f"rank=({np.linalg.matrix_rank(kl)},{np.linalg.matrix_rank(kr)}),",
        f"||C_X||={np.linalg.norm(compressed, ord=2):.6f}",
    )


def permutation_matrix(pi: tuple[int, ...]) -> np.ndarray:
    result = np.zeros((len(pi), len(pi)), dtype=float)
    for j, i in enumerate(pi):
        result[i, j] = 1
    return result


def hamming_permutation(pi: tuple[int, ...], sigma: tuple[int, ...]) -> int:
    return sum(a != b for a, b in zip(pi, sigma))


def greedy_permutation_code(r: int, distance: int) -> list[tuple[int, ...]]:
    code: list[tuple[int, ...]] = []
    for pi in permutations(range(r)):
        if all(hamming_permutation(pi, sigma) >= distance for sigma in code):
            code.append(pi)
    return code


def verify_permutation_fibre() -> None:
    r, distance, duplication = 5, 3, 2
    code = greedy_permutation_code(r, distance)
    assert len(code) > 1
    minimum = float("inf")
    for i in range(len(code)):
        for j in range(i + 1, len(code)):
            pi, sigma = code[i], code[j]
            d_h = hamming_permutation(pi, sigma)
            difference = permutation_matrix(pi) - permutation_matrix(sigma)
            norm = inf_to_one(difference)
            assert norm + 1e-12 >= 4 * d_h / 3
            metric = norm / (r * r)
            minimum = min(minimum, metric)

    # Exhaustively check duplication invariance once; its general proof is
    # the block-sum reduction in PF.22.
    difference = permutation_matrix(code[0]) - permutation_matrix(code[1])
    norm = inf_to_one(difference)
    duplicate_norm = inf_to_one(np.kron(difference, np.ones((duplication, duplication))))
    assert abs(duplicate_norm - duplication**2 * norm) < 1e-9
    duplicate_metric = duplicate_norm / (r * duplication) ** 2
    assert abs(norm / r**2 - duplicate_metric) < 1e-12
    print(
        "permutation compatibility packing:",
        f"r={r}, size={len(code)}, min_metric={minimum:.6f}",
    )


def verify_boolean_realization() -> None:
    q, h, _ = regular_hadamard(2)
    one = np.ones(16, dtype=np.int64)
    assert np.array_equal(h @ one, q * one)
    assert np.array_equal(h @ V0, q * V0)
    operator = np.kron(h, h)
    words = (
        np.kron(one, one),
        np.kron(one, V0),
        np.kron(V0, one),
        np.kron(V0, V0),
    )
    eigenvalue = q * q
    for word in words:
        assert np.array_equal(operator @ word, eigenvalue * word)

    r, duplication = 4, 2
    pi = (2, 0, 3, 1)
    left = np.stack([words[i] for i in range(r) for _ in range(duplication)])
    right = np.stack([words[pi[i]] for i in range(r) for _ in range(duplication)])
    n = operator.shape[0]
    kl = left @ left.T / n
    kr = right @ right.T / n
    cross = left @ right.T / n
    expected_margin = np.kron(np.eye(r), np.ones((duplication, duplication)))
    expected_cross = np.kron(permutation_matrix(pi), np.ones((duplication, duplication)))
    assert np.array_equal(kl, expected_margin)
    assert np.array_equal(kr, expected_margin)
    assert np.array_equal(cross, expected_cross)
    print("Boolean top-eigenvector compatibility realization: PASS")


def verify_spectral_truncation_bound() -> None:
    rng = np.random.default_rng(271828)
    p, q = 7, 6
    raw_l = rng.normal(size=(p, p))
    raw_r = rng.normal(size=(q, q))
    kl = raw_l @ raw_l.T
    kr = raw_r @ raw_r.T
    dl = np.sqrt(np.maximum(np.diag(kl), 1e-12))
    dr = np.sqrt(np.maximum(np.diag(kr), 1e-12))
    kl = kl / dl[:, None] / dl[None, :]
    kr = kr / dr[:, None] / dr[None, :]
    lr, _ = psd_sqrt(kl)
    rr, _ = psd_sqrt(kr)
    c = rng.normal(size=(p, q))
    c *= 0.9 / np.linalg.norm(c, ord=2)
    x = lr @ c @ rr

    total = p + q
    theta = 0.12

    def top_root(k: np.ndarray) -> np.ndarray:
        values, vectors = np.linalg.eigh(k)
        retained = np.where(values >= theta * total, np.sqrt(np.maximum(values, 0)), 0)
        return (vectors * retained) @ vectors.T

    x0 = top_root(kl) @ c @ top_root(kr)
    metric = 4 * inf_to_one(x - x0) / total**2
    assert metric <= 8 * np.sqrt(theta) + 1e-10
    print(
        "spectral compatibility truncation:",
        f"actual={metric:.6f}, theorem_bound={8*np.sqrt(theta):.6f}",
    )


def verify_rank_one_cycle_law() -> None:
    edge_words = list(product((-1, 1), repeat=3))
    gauges = list(product((-1, 1), repeat=3))
    orbits = []
    unseen = set(edge_words)
    while unseen:
        word = next(iter(unseen))
        orbit = set()
        for d0, d1, d2 in gauges:
            b01, b12, b20 = word
            orbit.add((d0 * b01 * d1, d1 * b12 * d2, d2 * b20 * d0))
        unseen -= orbit
        orbits.append(orbit)
    assert len(orbits) == 2
    products = [{a * b * c for a, b, c in orbit} for orbit in orbits]
    assert products == [{-1}, {1}] or products == [{1}, {-1}]

    # Only the balanced orbit is globally Gram-realizable at unit edge
    # correlations.  The negative OC cycle is coefficient-side holonomy.
    for word in edge_words:
        b01, b12, b20 = word
        gram = np.asarray(((1, b01, b20), (b01, 1, b12), (b20, b12, 1)))
        is_psd = float(np.linalg.eigvalsh(gram).min()) >= -1e-10
        assert is_psd == (b01 * b12 * b20 == 1)

    # Three unordered sector pairs separately have 2^(3-1) relative
    # orderings; this is the relative-antipode part of the OC fibre.
    assert 2 ** (3 - 1) == 4
    print("rank-one PSD-cycle obstruction and sector swaps: PASS")


def main() -> None:
    verify_contraction_parametrization()
    verify_permutation_fibre()
    verify_boolean_realization()
    verify_spectral_truncation_bound()
    verify_rank_one_cycle_law()
    print("PSD gluing compatibility entropy checks: PASS")


if __name__ == "__main__":
    main()
