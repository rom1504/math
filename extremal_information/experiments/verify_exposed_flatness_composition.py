#!/usr/bin/env python3
"""Finite checks for the exposed-flatness composition law."""

from __future__ import annotations

from math import sqrt

import numpy as np


def flatness(x: np.ndarray) -> float:
    return 1 - np.sum(np.abs(x)) / len(x)


def one_level_terms(blocks: list[np.ndarray]):
    N = sum(len(x) for x in blocks)
    u = np.concatenate(blocks)
    assert abs(float(u @ u) - N) < 1e-8
    allocation = 0.0
    transported = 0.0
    for x in blocks:
        lam = len(x) / N
        rho = np.linalg.norm(x) / sqrt(len(x))
        allocation += lam * (rho - 1) ** 2 / 2
        if rho > 0:
            transported += lam * rho * flatness(x / rho)
    return flatness(u), allocation, transported


def random_normalized_blocks(sizes, rng):
    blocks = [rng.normal(size=n) for n in sizes]
    norm = sqrt(sum(float(x @ x) for x in blocks))
    scale = sqrt(sum(sizes)) / norm
    return [scale * x for x in blocks]


def verify_one_level() -> int:
    rng = np.random.default_rng(260817)
    checks = 0
    for sizes in ((1, 3), (2, 5, 7), (4, 4, 4, 4)):
        for _ in range(100):
            blocks = random_normalized_blocks(sizes, rng)
            total, allocation, transported = one_level_terms(blocks)
            assert abs(total - allocation - transported) < 1e-10
            checks += 1
    return checks


def tree_expansion(u: np.ndarray, depth: int):
    N = len(u)
    root_rms = np.linalg.norm(u) / sqrt(N)
    internal_sum = 0.0
    current = [(u, 1.0)]
    level_weight_sums = []
    for _ in range(depth):
        nxt = []
        weight_sum = 0.0
        for block, global_weight in current:
            mid = len(block) // 2
            children = (block[:mid], block[mid:])
            parent_rms = np.linalg.norm(block) / sqrt(len(block))
            local_a = 0.0
            for child in children:
                lam = len(child) / len(block)
                rho = (
                    np.linalg.norm(child) / sqrt(len(child)) / parent_rms
                    if parent_rms > 0
                    else 0.0
                )
                local_a += lam * (rho - 1) ** 2 / 2
                child_global_rms = np.linalg.norm(child) / sqrt(len(child)) / root_rms
                child_weight = len(child) / N * child_global_rms
                nxt.append((child, child_weight))
            internal_sum += global_weight * local_a
            weight_sum += global_weight
        level_weight_sums.append(weight_sum)
        current = nxt
    leaf_sum = 0.0
    for leaf, weight in current:
        rms = np.linalg.norm(leaf) / sqrt(len(leaf))
        if rms > 0:
            leaf_sum += weight * flatness(leaf / rms)
    return internal_sum + leaf_sum, level_weight_sums


def verify_tree() -> int:
    rng = np.random.default_rng(817260)
    checks = 0
    for depth in range(1, 6):
        N = 2**depth
        for _ in range(50):
            u = rng.normal(size=N)
            u *= sqrt(N) / np.linalg.norm(u)
            expanded, weights = tree_expansion(u, depth)
            assert abs(expanded - flatness(u)) < 1e-10
            assert all(weight <= 1 + 1e-10 for weight in weights)
            checks += 1
    return checks


def verify_recovery() -> int:
    rng = np.random.default_rng(172608)
    checks = 0
    for N in (8, 16, 32):
        for _ in range(100):
            Q, _ = np.linalg.qr(rng.normal(size=(N, N)))
            eig = rng.uniform(-2, 2, size=N)
            M = Q @ np.diag(eig) @ Q.T
            Lambda = np.linalg.norm(M, 2)
            h = rng.normal(size=N)
            u = rng.normal(size=N)
            u *= sqrt(N) / np.linalg.norm(u)
            x = np.where(u >= 0, 1.0, -1.0)
            sphere = float(u @ M @ u) / 2 + float(h @ u)
            cube = float(x @ M @ x) / 2 + float(h @ x)
            kappa = np.linalg.norm(h) / (Lambda * sqrt(N))
            bound = Lambda * N * (1 + kappa) * sqrt(2 * flatness(u))
            assert sphere - cube <= bound + 1e-8
            checks += 1
    return checks


def tensor_amplitudes(delta: float, depth: int) -> np.ndarray:
    rho = (sqrt(1 + delta), sqrt(1 - delta))
    u = np.asarray([1.0])
    for _ in range(depth):
        u = np.kron(u, rho)
    return u


def verify_pumpable_tensor() -> int:
    checks = 0
    for delta in (0.1, 0.4, 0.8):
        s = (sqrt(1 + delta) + sqrt(1 - delta)) / 2
        for depth in range(1, 13):
            u = tensor_amplitudes(delta, depth)
            N = 2**depth
            assert abs(float(u @ u) - N) < 1e-9
            assert abs(np.sum(np.abs(u)) / N - s**depth) < 1e-12
            assert abs(flatness(u) - (1 - s**depth)) < 1e-12
            # Pure linear landscape: sphere N, cube ||u||_1.
            sphere = sqrt(N) * np.linalg.norm(u)
            cube = np.sum(np.abs(u))
            assert abs((sphere - cube) / N - (1 - s**depth)) < 1e-12
            expanded, weights = tree_expansion(u, depth)
            assert abs(expanded - (1 - s**depth)) < 1e-10
            assert all(weight <= 1 + 1e-10 for weight in weights)
            checks += 6
    return checks


def main() -> None:
    checks = verify_one_level()
    checks += verify_tree()
    checks += verify_recovery()
    checks += verify_pumpable_tensor()
    print(f"exposed-flatness composition checks passed: {checks}")


if __name__ == "__main__":
    main()
