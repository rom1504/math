#!/usr/bin/env python3
"""Exact finite checks for drafts/bcx_two_port_holonomy.md."""

from __future__ import annotations

from itertools import combinations, product
from math import sqrt
import random

import numpy as np


def parity(x: int) -> int:
    return bin(x).count("1") & 1


def walsh(d: int) -> np.ndarray:
    n = 1 << d
    return np.asarray(
        [[1 if parity(i & j) == 0 else -1 for j in range(n)] for i in range(n)],
        dtype=np.int64,
    )


def bent(m: int) -> np.ndarray:
    q = 1 << m
    return np.asarray(
        [1 if parity(u & v) == 0 else -1 for u in range(q) for v in range(q)],
        dtype=np.int64,
    )


def regular_hadamard(m: int) -> tuple[int, np.ndarray, np.ndarray]:
    q = 1 << m
    w = walsh(2 * m)
    b = bent(m)
    h = b[:, None] * w * b[None, :]
    a = h.copy()
    np.fill_diagonal(a, 0)
    n = q * q
    assert np.array_equal(h @ h, n * np.eye(n, dtype=np.int64))
    assert np.array_equal(h @ np.ones(n, dtype=np.int64), q * np.ones(n, dtype=np.int64))
    assert int(np.trace(h)) == 0
    return q, h, a


def cube(n: int) -> np.ndarray:
    return np.asarray(list(product((-1, 1), repeat=n)), dtype=np.int64)


def energies(a: np.ndarray, xs: np.ndarray) -> np.ndarray:
    return np.einsum("bi,ij,bj->b", xs, a, xs, optimize=True) // 2


def one_port_cap(
    old_energies: np.ndarray, xs: np.ndarray, t: np.ndarray, q: int
) -> int:
    best = 0
    fields = xs @ t
    for h, field in zip(old_energies, fields):
        for total in range(-q, q + 1, 2):
            clique = (total * total - q) // 2
            best = max(best, abs(int(h + field * total + clique)))
    return best


def check_one_port(m: int) -> None:
    q, h, a = regular_hadamard(m)
    n = q * q
    xs = cube(n)
    old = energies(a, xs)
    # Exhaust all queries at n=4.  At n=16 use every MM query, which is the
    # explicit query source used by the BCX regression.
    if n == 4:
        queries = xs
    else:
        queries = []
        for mask in range(1 << q):
            queries.append(
                np.asarray(
                    [
                        1 if (parity(u & v) ^ ((mask >> v) & 1)) == 0 else -1
                        for u in range(q)
                        for v in range(q)
                    ],
                    dtype=np.int64,
                )
            )
    max_gap = 0
    for t in queries:
        fp = one_port_cap(old, xs, t, q)
        fm = one_port_cap(-old, xs, t, q)
        gap = abs(fp - fm)
        assert gap <= q * (q - 1)
        max_gap = max(max_gap, gap)
    print(f"m={m}, n={n}: one-port queries={len(queries)}, max orientation gap={max_gap}")


def check_two_port(m: int) -> None:
    q, h, a = regular_hadamard(m)
    n = q * q
    tplus = np.asarray([[1, 1], [1, 1]], dtype=np.int64)
    tminus = np.asarray([[1, 1], [1, -1]], dtype=np.int64)
    mplus = np.kron(tplus, h)
    mminus = np.kron(tminus, h)
    assert np.array_equal(mminus @ mminus, 2 * n * np.eye(2 * n, dtype=np.int64))
    assert abs(float(np.linalg.norm(mplus.astype(float), ord=2)) - 2 * q) < 1e-8
    assert abs(float(np.linalg.norm(mminus.astype(float), ord=2)) - sqrt(2) * q) < 1e-8

    one = np.ones(n, dtype=np.int64)
    planted = int((one @ h @ one) // 2 + (one @ h @ one) // 2 + one @ h @ one)
    assert planted == 2 * q * n

    # Exhaust the complete two-shore Boolean cube at n=4.
    if n == 4:
        zs = cube(2 * n)
        vals_plus = np.einsum("bi,ij,bj->b", zs, mplus, zs, optimize=True) // 2
        vals_minus = np.einsum("bi,ij,bj->b", zs, mminus, zs, optimize=True) // 2
        qp = int(np.max(np.abs(vals_plus)))
        qm = int(np.max(np.abs(vals_minus)))
        assert qp == 2 * q * n
        assert qm <= sqrt(2) * q * n + 1e-8
        assert qp - qm >= (2 - sqrt(2)) * q * n - 1e-8
        print(f"m={m}, n={n}: exact joint caps Q+={qp}, Q-={qm}")
    else:
        print(
            f"m={m}, n={n}: spectral joint bounds Q+={2*q*n}, "
            f"Q-<={sqrt(2)*q*n:.6f}"
        )


def check_narrow_congruence() -> None:
    rng = random.Random(20260817)
    n, m = 4, 3
    xs, ys = cube(n), cube(m)
    for _ in range(24):
        a = np.zeros((n, n), dtype=np.int64)
        c = np.zeros((m, m), dtype=np.int64)
        for i, j in combinations(range(n), 2):
            a[i, j] = a[j, i] = rng.choice((-1, 1))
        for i, j in combinations(range(m), 2):
            c[i, j] = c[j, i] = rng.choice((-1, 1))
        b = np.asarray(
            [[rng.choice((-1, 1)) for _ in range(m)] for _ in range(n)],
            dtype=np.int64,
        )
        ha, hc = energies(a, xs), energies(c, ys)
        cross = xs @ b @ ys.T
        rp = int(np.max(np.abs(ha[:, None] + cross + hc[None, :])))
        rm = int(np.max(np.abs(-ha[:, None] + cross + hc[None, :])))
        qc = int(np.max(np.abs(hc)))
        assert abs(rp - rm) <= 2 * qc
    print("narrow-future orientation congruence: 24 exact instances passed")


def main() -> None:
    for m in (1, 2):
        check_one_port(m)
        check_two_port(m)
    check_narrow_congruence()
    print("BCX two-port holonomy checks passed")


if __name__ == "__main__":
    main()
