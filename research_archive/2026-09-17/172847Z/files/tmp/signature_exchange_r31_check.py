#!/usr/bin/env python3
"""Exact/numerical checks for the Wave 31 signature forest memo."""

from itertools import product

import numpy as np


def energy(a, x):
    return float(x @ a @ x)


def row_cost(a, x):
    y = a @ x
    return float(y @ y)


def signatures(words):
    base = words[0]
    return {tuple(int(words[j][i] * base[i]) for j in range(1, len(words)))
            for i in range(len(base))}


def projective_distance(x, y):
    d = int(np.count_nonzero(x != y))
    return min(d, len(x) - d)


def check_weighted_separation():
    a = np.zeros((5, 5), dtype=float)
    for i, j, value in [
        (0, 1, 10), (0, 2, 10), (1, 2, 20), (3, 4, 20),
        (0, 3, 1), (0, 4, 1), (1, 3, -1), (1, 4, -1),
    ]:
        a[i, j] = a[j, i] = value

    states = [np.array(x, dtype=float) for x in product((-1, 1), repeat=5)]
    values = [energy(a, x) for x in states]
    q = max(values)
    x = np.ones(5)
    xf = x.copy()
    xf[:3] *= -1
    assert q == 120.0
    assert energy(a, x) == q == energy(a, xf)
    assert row_cost(a, x) == 2968.0
    assert row_cost(a, xf) == 3048.0

    w = a * np.outer(x, x)
    shore = np.array([1, 1, 1, 0, 0], dtype=bool)
    avec = np.array([
        sum(w[i, j] for j in range(5) if shore[i] != shore[j])
        for i in range(5)
    ])
    rvec = w.sum(axis=1)
    assert np.array_equal(avec, np.array([2, -2, 0, 0, 0]))
    rhs = row_cost(a, x) + 4 * np.sum(avec * avec - rvec * avec)
    assert rhs == row_cost(a, xf)


def check_random_forest(seed=3101, trials=100):
    rng = np.random.default_rng(seed)
    for _ in range(trials):
        n = int(rng.integers(5, 11))
        upper = rng.choice((-1.0, 1.0), size=(n, n))
        a = np.triu(upper, 1)
        a += a.T
        root = rng.choice((-1.0, 1.0), size=n)

        # A star is enough to check the theorem.  Orient every child so the
        # displayed Hamming distance is projectively minimal.
        q = int(rng.integers(1, 6))
        words = [root]
        dsum = 0
        union = set()
        for _j in range(q):
            y = rng.choice((-1.0, 1.0), size=n)
            if np.count_nonzero(y != root) > n / 2:
                y = -y
            diff = set(np.flatnonzero(y != root).tolist())
            dsum += len(diff)
            union |= diff
            words.append(y)

        assert len(signatures(words)) <= 1 + len(union) <= 1 + dsum

        # Enumerate the coarsest signature coset exactly for these small n.
        sig_to_indices = {}
        base = words[0]
        for i in range(n):
            sig = tuple(int(words[j][i] * base[i]) for j in range(1, len(words)))
            sig_to_indices.setdefault(sig, []).append(i)
        classes = list(sig_to_indices.values())
        c_sig = 0.0
        for signs in product((-1.0, 1.0), repeat=len(classes)):
            h = base.copy()
            for sign, inds in zip(signs, classes):
                h[inds] = base[inds] * sign
            c_sig = max(c_sig, row_cost(a, h))

        op2 = float(np.linalg.norm(a, 2) ** 2)
        bound = 2 * row_cost(a, root) + 8 * op2 * dsum
        assert c_sig <= bound + 1e-7


if __name__ == "__main__":
    check_weighted_separation()
    check_random_forest()
    print("PASS: exchange susceptibility, weighted separation, and forest closure")
