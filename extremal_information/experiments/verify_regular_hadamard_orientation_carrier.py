#!/usr/bin/env python3
"""Exact finite checks for the regular-Hadamard orientation carrier."""

from __future__ import annotations

import itertools
from collections import defaultdict

import numpy as np


def signs(k: int):
    return itertools.product((-1, 1), repeat=k)


def edges(k: int):
    return list(itertools.combinations(range(k), 2))


def projective_switch_orbit(T: np.ndarray) -> set[tuple[int, ...]]:
    k = len(T)
    out: set[tuple[int, ...]] = set()
    for d0 in signs(k):
        d = np.array(d0, dtype=int)
        for epsilon in (-1, 1):
            U = epsilon * d[:, None] * T * d[None, :]
            out.add(tuple(int(v) for v in U.ravel()))
    return out


def scalar_matrix(sigmas, bridge_signs, support=None):
    k = len(sigmas)
    ee = edges(k)
    if support is None:
        support = (1,) * len(ee)
    T = np.diag(np.array(sigmas, dtype=int))
    for a, b, (i, j) in zip(support, bridge_signs, ee):
        if a:
            T[i, j] = T[j, i] = b
    return T


def exact_cap(T: np.ndarray, H: np.ndarray) -> int:
    N = len(T) * len(H)
    M = np.kron(T, H)
    best = 0
    for z0 in signs(N):
        z = np.array(z0, dtype=int)
        val = abs(int(z @ M @ z)) // 2
        best = max(best, val)
    return best


def charpoly_coefficients(T: np.ndarray) -> tuple[int, ...]:
    """Exact Faddeev--LeVerrier coefficients for a small integer matrix."""
    k = len(T)
    A = T.astype(object)
    B = np.eye(k, dtype=object)
    coeffs = [1]
    for r in range(1, k + 1):
        B = A @ B
        c = -sum(B[i, i] for i in range(k)) // r
        coeffs.append(int(c))
        B = B + c * np.eye(k, dtype=object)
    return tuple(coeffs)


def connected_components(k: int, support) -> int:
    ee = edges(k)
    adj = [[] for _ in range(k)]
    for keep, (i, j) in zip(support, ee):
        if keep:
            adj[i].append(j)
            adj[j].append(i)
    seen = set()
    count = 0
    for root in range(k):
        if root in seen:
            continue
        count += 1
        stack = [root]
        seen.add(root)
        while stack:
            i = stack.pop()
            for j in adj[i]:
                if j not in seen:
                    seen.add(j)
                    stack.append(j)
    return count


def check_orbit_counts() -> int:
    checks = 0
    for k in range(1, 5):
        ee = edges(k)
        for support in itertools.product((0, 1), repeat=len(ee)):
            e = sum(support)
            c = connected_components(k, support)
            representatives = {}
            for sigma in signs(k):
                active = [a for a, keep in enumerate(support) if keep]
                for active_signs in signs(len(active)):
                    bridge = [1] * len(ee)
                    for a, b in zip(active, active_signs):
                        bridge[a] = b
                    T = scalar_matrix(sigma, bridge, support)
                    key = min(projective_switch_orbit(T))
                    representatives[key] = True
            expected = 2 ** (e + c - 1)
            assert len(representatives) == expected, (k, support, expected)
            checks += 1
    return checks


def check_gauge_invariance(H: np.ndarray) -> int:
    checks = 0
    rng = np.random.default_rng(260817)
    for k in (2, 3):
        for _ in range(12):
            sigma = tuple(int(x) for x in rng.choice((-1, 1), size=k))
            bridge = tuple(int(x) for x in rng.choice((-1, 1), size=len(edges(k))))
            T = scalar_matrix(sigma, bridge)
            cap = exact_cap(T, H)
            d = rng.choice((-1, 1), size=k)
            for epsilon in (-1, 1):
                U = epsilon * d[:, None] * T * d[None, :]
                assert exact_cap(U, H) == cap
                checks += 1
    return checks


def check_gluing_fibres() -> int:
    """Count fibres for singleton pieces: r-edge connected simple graphs."""
    checks = 0
    # Singleton marginal projective carriers are unique.  For every connected
    # graph the joined carrier count must be 2^r.
    for k in range(2, 5):
        ee = edges(k)
        for support in itertools.product((0, 1), repeat=len(ee)):
            if connected_components(k, support) != 1:
                continue
            r = sum(support)
            carriers = set()
            active = [a for a, keep in enumerate(support) if keep]
            for sigma in signs(k):
                for active_signs in signs(len(active)):
                    bridge = [1] * len(ee)
                    for a, b in zip(active, active_signs):
                        bridge[a] = b
                    T = scalar_matrix(sigma, bridge, support)
                    carriers.add(min(projective_switch_orbit(T)))
            assert len(carriers) == 2**r, (k, support, r, len(carriers))
            checks += 1
    return checks


def main() -> None:
    H = np.array(
        [
            [1, 1, 1, -1],
            [1, -1, 1, 1],
            [1, 1, -1, 1],
            [-1, 1, 1, 1],
        ],
        dtype=int,
    )
    assert np.array_equal(H @ H, 4 * np.eye(4, dtype=int))
    assert np.array_equal(H @ np.ones(4, dtype=int), 2 * np.ones(4, dtype=int))
    assert np.trace(H) == 0

    checks = check_orbit_counts()
    checks += check_gauge_invariance(H)
    checks += check_gluing_fibres()

    T_plus = np.array([[1, 1], [1, 1]], dtype=int)
    T_minus = np.array([[1, 1], [1, -1]], dtype=int)
    assert np.array_equal(T_minus @ T_minus, 2 * np.eye(2, dtype=int))
    assert exact_cap(T_plus, H) == 16
    assert exact_cap(T_minus, H) == 10
    checks += 4

    T_balanced = np.ones((3, 3), dtype=int)
    T_unbalanced = np.array([[1, 1, 1], [1, 1, -1], [1, -1, 1]], dtype=int)
    assert charpoly_coefficients(T_unbalanced) == (1, -3, 0, 4)
    assert np.array_equal(
        T_unbalanced @ T_unbalanced @ T_unbalanced
        - 3 * T_unbalanced @ T_unbalanced
        + 4 * np.eye(3, dtype=int),
        np.zeros((3, 3), dtype=int),
    )
    assert exact_cap(T_balanced, H) == 36
    assert exact_cap(T_unbalanced, H) == 20
    checks += 4

    T0 = np.array(
        [
            [-1, -1, 0, -1],
            [-1, -1, 0, -1],
            [0, 0, -1, 0],
            [-1, -1, 0, 1],
        ],
        dtype=int,
    )
    T1 = np.array(
        [
            [-1, -1, -1, -1],
            [-1, -1, 0, 0],
            [-1, 0, -1, 0],
            [-1, 0, 0, 1],
        ],
        dtype=int,
    )
    expected_poly = (1, 2, -3, -4, 0)
    assert charpoly_coefficients(T0) == expected_poly
    assert charpoly_coefficients(T1) == expected_poly
    assert exact_cap(T0, H) == 32
    assert exact_cap(T1, H) == 34
    assert not (projective_switch_orbit(T0) & projective_switch_orbit(T1))
    checks += 5

    print(f"regular-Hadamard orientation-carrier checks passed: {checks}")


if __name__ == "__main__":
    main()
