#!/usr/bin/env python3
"""Finite checks for common-pole synchronization and its tensor law."""

from __future__ import annotations

from itertools import product

import numpy as np

from verify_bcx_two_port_holonomy import regular_hadamard


def cube(n: int) -> np.ndarray:
    return np.asarray(list(product((-1, 1), repeat=n)), dtype=np.int64)


def exact_boolean(H: np.ndarray, ports: np.ndarray, m: int) -> int:
    xs = cube(len(H))
    child = np.abs(np.einsum("bi,ij,bj->b", xs, H, xs) // 2)
    fields = m * np.sum(np.abs(xs @ ports.T), axis=1)
    return int(np.max(child + fields))


def spherical_upper(r: int, n: int, p: int, m: int) -> int:
    return r * n // 2 + m * p * n


def deficit(pole: np.ndarray, ports: np.ndarray) -> float:
    p, n = ports.shape
    return 1.0 - float(np.sum(np.abs(ports @ pole))) / (p * n)


def finite_recovery_checks() -> int:
    rng = np.random.default_rng(817261)
    checks = 0
    q, H, _ = regular_hadamard(1)
    n = len(H)
    pole = np.ones(n, dtype=np.int64)
    xs = cube(n)
    top = xs[(xs @ H.T == q * xs).all(axis=1)]
    assert len(top) > 0
    for p in (1, 2, 3):
        for m in (1, 2):
            for _ in range(20):
                ports = rng.choice((-1, 1), size=(p, n)).astype(np.int64)
                B = exact_boolean(H, ports, m)
                S_upper = spherical_upper(q, n, p, m)
                delta = deficit(pole, ports)
                bound = (m * p / q) * delta * q * n
                assert S_upper - B <= bound + 1e-9
                checks += 1
    return checks


def tensor_checks() -> int:
    rng = np.random.default_rng(817262)
    checks = 0
    for n1, n2, p1, p2 in ((4, 4, 3, 2), (4, 8, 2, 3), (8, 8, 4, 3)):
        x1 = rng.choice((-1, 1), size=n1).astype(np.int64)
        x2 = rng.choice((-1, 1), size=n2).astype(np.int64)
        W1 = rng.choice((-1, 1), size=(p1, n1)).astype(np.int64)
        W2 = rng.choice((-1, 1), size=(p2, n2)).astype(np.int64)
        d1 = deficit(x1, W1)
        d2 = deficit(x2, W2)
        W = np.asarray([np.kron(w, v) for w in W1 for v in W2])
        d = deficit(np.kron(x1, x2), W)
        assert abs((1 - d) - (1 - d1) * (1 - d2)) < 1e-12
        assert d <= d1 + d2 + 1e-12
        checks += 2
    return checks


def completion_lipschitz_checks() -> int:
    # Direct generic max-Lipschitz check on a tiny joined problem.
    rng = np.random.default_rng(817263)
    checks = 0
    for old_n, aux_n in ((3, 2), (4, 3)):
        all_x = cube(old_n + aux_n)
        for _ in range(20):
            old = rng.integers(-4, 5, size=len(all_x))
            aux = rng.integers(-3, 4, size=len(all_x))
            base = int(np.max(old))
            joined = int(np.max(old + aux))
            assert abs(joined - base) <= int(np.max(np.abs(aux)))
            checks += 1
    return checks


def main() -> None:
    checks = finite_recovery_checks()
    checks += tensor_checks()
    checks += completion_lipschitz_checks()
    print(f"common-pole synchronization checks passed: {checks}")


if __name__ == "__main__":
    main()
