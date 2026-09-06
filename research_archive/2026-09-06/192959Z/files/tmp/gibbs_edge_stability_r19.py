#!/usr/bin/env python3
"""Exact finite-state checks for Gibbs edge/star variation identities."""

from __future__ import annotations

import itertools

import numpy as np


def states(n: int):
    for sigma in (-1, 1):
        for x in itertools.product((-1, 1), repeat=n):
            yield sigma, np.array(x, dtype=float)


def z(A: np.ndarray, beta: float) -> float:
    return sum(np.exp(beta * sigma * (x @ A @ x)) for sigma, x in states(len(A)))


def audit(A: np.ndarray, beta: float, i: int) -> None:
    n = len(A)
    base = z(A, beta)

    # Edge identity.
    j = (i + 1) % n
    Ae = A.copy()
    Ae[i, j] *= -1
    Ae[j, i] *= -1
    moment = 0.0
    for sigma, x in states(n):
        weight = np.exp(beta * sigma * (x @ A @ x)) / base
        moment += weight * np.exp(-4 * beta * sigma * A[i, j] * x[i] * x[j])
    assert np.isclose(z(Ae, beta) / base, moment)

    # Row and biased-star identities.
    keep = [j for j in range(n) if j != i]
    C = A[np.ix_(keep, keep)]
    b = A[i, keep]
    child_z = z(C, beta)

    def F(c: np.ndarray) -> float:
        value = 0.0
        for sigma, y in states(n - 1):
            weight = np.exp(beta * sigma * (y @ C @ y)) / child_z
            value += weight * np.cosh(2 * beta * (c @ y))
        return value

    assert np.isclose(base / (2 * child_z), F(b))
    rows = [np.array(c, dtype=float) for c in itertools.product((-1, 1), repeat=n - 1)]
    assert np.isclose(np.mean([F(c) for c in rows]), np.cosh(2 * beta) ** (n - 1))

    for p in (0.2, 0.5, 0.8):
        lhs = 0.0
        for bits in itertools.product((0, 1), repeat=n - 1):
            flips = np.array(bits, dtype=bool)
            prob = p ** sum(bits) * (1 - p) ** (n - 1 - sum(bits))
            c = b.copy()
            c[flips] *= -1
            lhs += prob * F(c) / F(b)
        ap = np.sqrt(1 + 2 * p * (1 - p) * (np.cosh(4 * beta) - 1))
        lam = 0.5 * np.log((1 - p + p * np.exp(4 * beta)) / (1 - p + p * np.exp(-4 * beta)))
        rhs = ap ** (n - 1) * F((1 - lam / (2 * beta)) * b) / F(b)
        assert np.isclose(lhs, rhs)


def main() -> None:
    A = np.array([[0, 1, -1, 1], [1, 0, 1, -1], [-1, 1, 0, 1], [1, -1, 1, 0]], dtype=float)
    for beta in (0.07, 0.3, 1.1):
        for i in range(len(A)):
            audit(A, beta, i)
    print("PASS: edge, row, uniform-star, and Hamming-p transforms")


if __name__ == "__main__":
    main()
