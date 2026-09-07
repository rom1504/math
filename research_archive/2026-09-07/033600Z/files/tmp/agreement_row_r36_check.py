#!/usr/bin/env python3
"""Finite checks for the Wave 36 anchored row-transfer memo."""

from itertools import combinations, product
import math
import random

import numpy as np


def signing(n: int, seed: int) -> np.ndarray:
    rng = random.Random(seed)
    a = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(i + 1, n):
            a[i, j] = a[j, i] = rng.choice((-1, 1))
    return a


def row_square(a: np.ndarray, z: np.ndarray) -> int:
    return int(np.dot(a @ z, a @ z))


def local_row_square(a: np.ndarray, z: np.ndarray, s: tuple[int, ...]) -> int:
    idx = np.array(s, dtype=int)
    b = a[np.ix_(idx, idx)]
    y = z[idx]
    return int(np.dot(b @ y, b @ y))


def q_value(a: np.ndarray) -> tuple[int, np.ndarray, int]:
    best = -1
    best_x = None
    best_sigma = None
    for bits in product((-1, 1), repeat=len(a)):
        x = np.array(bits, dtype=np.int64)
        e = int(x @ a @ x)
        if abs(e) > best:
            best = abs(e)
            best_x = x
            best_sigma = 1 if e >= 0 else -1
    assert best_x is not None and best_sigma is not None
    return best, best_x, best_sigma


def check_uniform_and_anchored_identity() -> None:
    n, m, v = 8, 5, 2
    a = signing(n, 3601)
    z = np.array([1, -1, 1, 1, -1, -1, 1, -1], dtype=np.int64)
    subsets = list(combinations(range(n), m))
    actual = sum(local_row_square(a, z, s) for s in subsets) / len(subsets)
    p2 = m * (m - 1) / (n * (n - 1))
    p3 = m * (m - 1) * (m - 2) / (n * (n - 1) * (n - 2))
    predicted = p3 * row_square(a, z) + (p2 - p3) * n * (n - 1)
    assert abs(actual - predicted) < 1e-10

    w = np.diag(z) @ a @ np.diag(z)
    r = w @ np.ones(n, dtype=np.int64)
    triples_v = []
    for j, k in combinations([i for i in range(n) if i != v], 2):
        triples_v.append(
            int(w[v, j] * w[v, k] + w[v, j] * w[j, k] + w[v, k] * w[j, k])
        )
    c_v = sum(triples_v)
    c_v_formula = (r[v] ** 2) / 2 + (w @ r)[v] - 3 * (n - 1) / 2
    assert c_v == c_v_formula

    anchored = [s for s in subsets if v in s]
    actual_v = sum(local_row_square(a, z, s) for s in anchored) / len(anchored)
    alpha2 = (m - 1) * (m - 2) / ((n - 1) * (n - 2))
    alpha3 = (m - 1) * (m - 2) * (m - 3) / ((n - 1) * (n - 2) * (n - 3))
    predicted_v = (
        m * (m - 1)
        + alpha3 * (row_square(a, z) - n * (n - 1))
        + 2 * (alpha2 - alpha3) * c_v
    )
    assert abs(actual_v - predicted_v) < 1e-10


def check_ground_stability_and_join() -> None:
    # Gauge a small arbitrary signing at an absolute ground.
    b0 = signing(6, 3602)
    qb, x, sigma = q_value(b0)
    b = sigma * np.diag(x) @ b0 @ np.diag(x)
    one_l = np.ones(len(b), dtype=np.int64)
    assert int(one_l @ b @ one_l) == qb

    rows = b @ one_l
    assert np.all(rows >= 0)
    assert int(rows @ rows) <= (len(b) - 1) * qb
    for bits in product((0, 1), repeat=len(b)):
        f = np.array(bits, dtype=bool)
        cut = int(b[np.ix_(f, ~f)].sum())
        assert 0 <= cut <= qb // 2

    # Adjoin positive hubs.  The old ground remains an exact absolute ground.
    h, ell = 2, len(b)
    n = h + ell
    a = np.ones((n, n), dtype=np.int64)
    np.fill_diagonal(a, 0)
    a[h:, h:] = b
    qa, xa, sigmaa = q_value(a)
    expected = qb + 2 * h * ell + h * (h - 1)
    assert qa == expected
    assert int(np.ones(n, dtype=np.int64) @ a @ np.ones(n, dtype=np.int64)) == qa
    assert row_square(a, np.ones(n, dtype=np.int64)) >= h * (n - 1) ** 2
    assert sigmaa in (-1, 1) and len(xa) == n


def check_hamming_bound() -> None:
    n, m = 9, 6
    a = signing(n, 3603)
    s = tuple(range(m))
    y = np.ones(n, dtype=np.int64)
    z = y.copy()
    z[[1, 4]] *= -1
    lhs = abs(math.sqrt(local_row_square(a, z, s)) - math.sqrt(local_row_square(a, y, s)))
    op = float(np.linalg.norm(a[np.ix_(s, s)], ord=2))
    rhs = 2 * op * math.sqrt(2)
    assert lhs <= rhs + 1e-10


if __name__ == "__main__":
    check_uniform_and_anchored_identity()
    check_ground_stability_and_join()
    check_hamming_bound()
    print("PASS: uniform/anchored identities, ground stability, hub join, and Hamming bound")
