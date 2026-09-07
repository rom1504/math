#!/usr/bin/env python3
"""Algebraic checker for the Wave-18 symplectic switching theorem."""

from __future__ import annotations

import itertools

import numpy as np


def dot2(x: int, y: int) -> int:
    return bin(x & y).count("1") & 1


def phase(bit: int) -> int:
    return -1 if bit else 1


def matrices(k: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    n = k * k
    K = np.empty((n, n), dtype=np.int64)
    for u, v, x, y in itertools.product(range(k), repeat=4):
        K[u * k + v, x * k + y] = phase(dot2(v, x) ^ dot2(u, y))
    P = np.zeros_like(K)
    for u in range(k):
        r = np.array([phase(dot2(u, v)) for v in range(k)], dtype=np.int64)
        sl = slice(u * k, (u + 1) * k)
        P[sl, sl] = np.outer(r, r)
    eye = np.eye(n, dtype=np.int64)
    A = K - eye
    D = P - eye
    C = K - P
    return K, P, A, C, D


def state(k: int, a: int, b: int) -> np.ndarray:
    return np.array(
        [phase(dot2(u, v) ^ dot2(a, u) ^ dot2(b, v)) for u in range(k) for v in range(k)],
        dtype=np.int64,
    )


def audit(k: int) -> None:
    n = k * k
    K, P, A, C, D = matrices(k)
    eye = np.eye(n, dtype=np.int64)
    assert np.array_equal(K, K.T)
    assert np.array_equal(np.diag(K), np.ones(n, dtype=np.int64))
    assert np.array_equal(K @ K, n * eye)
    assert np.array_equal(P @ P, k * P)
    assert np.array_equal(K @ P, k * P)
    assert np.array_equal(P @ K, k * P)
    assert np.array_equal(C @ C, k * k * eye - k * P)
    assert np.array_equal(A, C + D)

    b = 1
    a = 0
    c = 1
    assert dot2(a, b) == 0 and dot2(c, b) == 1
    z = state(k, a, b)
    w = state(k, a ^ c, b)
    g = np.repeat(np.array([phase(dot2(c, u)) for u in range(k)], dtype=np.int64), k)
    assert np.array_equal(w, z * g)
    assert np.array_equal(P @ z, np.zeros(n, dtype=np.int64))
    assert np.array_equal(P @ w, np.zeros(n, dtype=np.int64))
    assert np.array_equal(C @ z, k * z)
    assert np.array_equal(C @ w, -k * w)
    assert np.array_equal(D @ z, -z)
    assert np.array_equal(D @ w, -w)
    assert int(z @ w) == 0

    qc = n * k
    qa = n * (k + 1)
    assert int(z @ C @ z) == qc
    assert int(w @ C @ w) == -qc
    assert int(z @ A @ z) == n * (k - 1)
    assert int(w @ A @ w) == -qa
    assert int(z @ A @ w) == 0

    same = np.flatnonzero(g == 1)
    diff = np.flatnonzero(g == -1)
    assert len(same) == len(diff) == n // 2
    h_same = int(z[same] @ A[np.ix_(same, same)] @ z[same])
    h_diff = int(z[diff] @ A[np.ix_(diff, diff)] @ z[diff])
    l_same = int(np.abs(A[np.ix_(diff, same)] @ z[same]).sum())
    l_diff = int(np.abs(A[np.ix_(same, diff)] @ z[diff]).sum())
    assert h_same == h_diff == -n // 2
    assert l_same == l_diff == n * k // 2
    raw = []
    for h, ell in ((h_same, l_same), (h_diff, l_diff)):
        raw.extend((max(2 * ell + h - qa, 0), max(2 * ell - h - qa, 0)))
    assert raw == [0, 0, 0, 0]

    print(
        f"k={k} n={n}: cells={k}x{k}, S={n * (k - 1)}, "
        f"Q(C)={qc}, Q(A)={qa}, orbit-gap=0, raw={tuple(raw)}"
    )


def main() -> None:
    for k in (2, 4, 8, 16):
        audit(k)
    print("PASS: symplectic Reynolds, eigenspace, antipodal, and raw-shore identities")


if __name__ == "__main__":
    main()
