#!/usr/bin/env python3
"""Numerical audit for the three-fibre conference obstruction.

Build A_k = C_k tensor (J_3-2I_3) + I_k tensor (J_3-I_3), then run
greedy/randomized coordinate ascent for max |x^T A_k x|.  This is a
heuristic evaluator; all algebraic identities are checked exactly.
"""

from __future__ import annotations

import math
import numpy as np
from scipy.stats import norm


def legendre(x: int, p: int) -> int:
    x %= p
    if x == 0:
        return 0
    return 1 if pow(x, (p - 1) // 2, p) == 1 else -1


def paley_conference(p: int) -> np.ndarray:
    assert p % 4 == 1
    k = p + 1
    c = np.zeros((k, k), dtype=np.int64)
    c[0, 1:] = 1
    c[1:, 0] = 1
    for x in range(p):
        for y in range(x + 1, p):
            c[x + 1, y + 1] = c[y + 1, x + 1] = legendre(x - y, p)
    assert np.array_equal(c @ c, (k - 1) * np.eye(k, dtype=np.int64))
    return c


def fibre_matrix(c: np.ndarray) -> np.ndarray:
    k = len(c)
    r = np.ones((3, 3), dtype=np.int64) - 2 * np.eye(3, dtype=np.int64)
    d = np.ones((3, 3), dtype=np.int64) - np.eye(3, dtype=np.int64)
    a = np.kron(c, r) + np.kron(np.eye(k, dtype=np.int64), d)
    assert np.all(np.diag(a) == 0)
    assert np.all(np.abs(a - np.diag(np.diag(a))) == 1 - np.eye(3 * k, dtype=np.int64))
    return a


def coordinate_ascent(a: np.ndarray, rng: np.random.Generator, sweeps: int = 100) -> tuple[int, np.ndarray]:
    n = len(a)
    x = rng.choice(np.array([-1, 1], dtype=np.int64), size=n)
    h = a @ x
    value = int(x @ h)
    for _ in range(sweeps):
        changed = False
        for i in rng.permutation(n):
            delta = -4 * int(x[i]) * int(h[i])
            if delta > 0:
                old = int(x[i])
                x[i] = -old
                h -= 2 * old * a[:, i]
                value += delta
                changed = True
        if not changed:
            break
    assert value == int(x @ a @ x)
    return value, x


def audit(p: int, restarts: int, seed: int = 1) -> None:
    c = paley_conference(p)
    a = fibre_matrix(c)
    n = len(a)
    rng = np.random.default_rng(seed)
    best = -10**30
    for orient in (1, -1):
        aa = orient * a
        for _ in range(restarts):
            value, _ = coordinate_ascent(aa, rng)
            best = max(best, value)
    print(
        f"k={len(c):4d} n={n:4d} best={best:10d} "
        f"ratio={best/(n*math.sqrt(n-1)):.9f}"
    )


def dependent_rounding_audit(p: int, trials: int, seed: int = 1) -> None:
    c = paley_conference(p)
    a = fibre_matrix(c)
    n = len(a)
    bmat = a / math.sqrt(n - 1)
    rng = np.random.default_rng(seed)
    vals = []
    t = 0.8414699114
    for _ in range(trials):
        s = rng.choice(np.array([-1.0, 1.0]), size=n)
        g = bmat @ s
        f = np.sign(g + t * s)
        f[f == 0] = 1
        y = np.sign(
            -0.1225227631 * s
            - 0.1063660230 * g
            + 0.4194582166 * f
            + 0.8158276928 * (bmat @ f)
        )
        y[y == 0] = 1
        vals.append(float(y @ bmat @ y) / n)
    print(
        f"AMP k={len(c):4d} n={n:4d} "
        f"mean={np.mean(vals):.9f} max={np.max(vals):.9f} "
        f"absmean={np.mean(np.abs(vals)):.9f}"
    )


def onsager_corrected_audit(p: int, trials: int, seed: int = 1) -> None:
    c = paley_conference(p)
    a_mat = fibre_matrix(c)
    n = len(a_mat)
    bmat = a_mat / math.sqrt(n - 1)
    rng = np.random.default_rng(seed)
    vals = []
    t = 0.8414699114
    aa = 2 * norm.pdf(t)
    bb = 2 * norm.cdf(t) - 1
    ss = math.sqrt(1 - aa * aa - bb * bb)
    p0, q0, r0, d0 = (
        0.5859761744,
        0.6179560304,
        0.2396817825,
        0.4661704739,
    )
    for _ in range(trials):
        spin = rng.choice(np.array([-1.0, 1.0]), size=n)
        g = bmat @ spin
        f = np.sign(g + t * spin)
        f[f == 0] = 1
        residual = f - aa * g - bb * spin
        w = bmat @ residual
        y = np.sign(p0 * spin + q0 * g + (r0 / ss) * residual + (d0 / ss) * w)
        y[y == 0] = 1
        vals.append(float(y @ bmat @ y) / n)
    print(
        f"COR k={len(c):4d} n={n:4d} "
        f"mean={np.mean(vals):.9f} max={np.max(vals):.9f} "
        f"absmean={np.mean(np.abs(vals)):.9f}"
    )


if __name__ == "__main__":
    for prime in (13, 17, 29, 37, 41, 61, 73, 89, 97, 101, 109, 149, 193):
        audit(prime, restarts=100)
        dependent_rounding_audit(prime, trials=200)
        onsager_corrected_audit(prime, trials=200)
