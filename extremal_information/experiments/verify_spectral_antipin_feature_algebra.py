#!/usr/bin/env python3
"""Finite checks for the spectral anti-pin and finite-port Gram algebra."""

from __future__ import annotations

from itertools import product
from math import sqrt

import numpy as np

from verify_bounded_cap_contextual_metric_compiler import build, matvec


V0 = (
    -1, -1, -1, 1,
    -1, -1, 1, -1,
    1, -1, 1, 1,
    -1, 1, 1, 1,
)


def hollow_energy(H: list[list[int]], x: tuple[int, ...]) -> int:
    # The audited matrices have trace zero, so this equals x^T H x / 2.
    n = len(x)
    return sum(H[i][j] * x[i] * x[j]
               for i in range(n) for j in range(i + 1, n))


def port_cap(H: list[list[int]], ports: tuple[tuple[int, ...], ...], r: int) -> int:
    best = 0
    for x in product((-1, 1), repeat=len(H)):
        value = abs(hollow_energy(H, x))
        value += r * sum(abs(sum(a * b for a, b in zip(w, x))) for w in ports)
        best = max(best, value)
    return best


def gram(H: np.ndarray, ports: list[np.ndarray], r: int):
    n = len(ports[0])
    G = np.asarray([[int(a @ b) / n for b in ports] for a in ports])
    R = np.asarray([[int(a @ H @ b) / (r * n) for b in ports]
                    for a in ports])
    return G, R


def verify_resolvents() -> None:
    q, n, H_list, _ = build(2)
    H = np.asarray(H_list, dtype=float)
    eye = np.eye(n)
    for w_tuple in (tuple([1] * n), V0, tuple((-1) ** i for i in range(n))):
        w = np.asarray(w_tuple, dtype=float)
        rho = float(w @ H @ w) / (q * n)
        psi_plus = q / (2 * n) * float(w @ np.linalg.solve(2 * q * eye - H, w))
        psi_minus = q / (2 * n) * float(w @ np.linalg.solve(2 * q * eye + H, w))
        assert abs(psi_plus - (2 + rho) / 6) < 1e-10
        assert abs(psi_minus - (2 - rho) / 6) < 1e-10


def verify_cross_gram_counterexample() -> None:
    q, n, H_list, _ = build(2)
    H = np.asarray(H_list, dtype=np.int64)
    one = np.ones(n, dtype=np.int64)
    v = np.asarray(V0, dtype=np.int64)
    assert int(one @ v) == 0
    assert np.array_equal(H @ one, q * one)
    assert np.array_equal(H @ v, q * v)

    G_same, R_same = gram(H, [one, one], q)
    G_cross, R_cross = gram(H, [one, v], q)
    assert np.array_equal(G_same, np.ones((2, 2)))
    assert np.array_equal(R_same, np.ones((2, 2)))
    assert np.array_equal(G_cross, np.eye(2))
    assert np.array_equal(R_cross, np.eye(2))

    one_tuple = tuple(int(x) for x in one)
    assert port_cap(H_list, (one_tuple,), q) == 3 * q * n // 2
    assert port_cap(H_list, (V0,), q) == 3 * q * n // 2
    same = port_cap(H_list, (one_tuple, one_tuple), q)
    cross = port_cap(H_list, (one_tuple, V0), q)
    assert same == 5 * q * n // 2
    assert cross <= (0.5 + sqrt(2)) * q * n + 1e-9
    assert same - cross >= (2 - sqrt(2)) * q * n - 1e-9
    print(
        f"n={n}: one-port={3*q*n//2}, same-pair={same}, "
        f"orthogonal-pair={cross}, gap={same-cross}"
    )

    # The tensor witness remains balanced and a top Boolean eigenvector.
    H2 = np.kron(H, H)
    one2 = np.ones(n * n, dtype=np.int64)
    v2 = np.kron(v, one)
    assert int(one2 @ v2) == 0
    assert np.array_equal(H2 @ one2, q * q * one2)
    assert np.array_equal(H2 @ v2, q * q * v2)


def main() -> None:
    verify_resolvents()
    verify_cross_gram_counterexample()
    print("spectral anti-pin feature-algebra checks: PASS")


if __name__ == "__main__":
    main()
