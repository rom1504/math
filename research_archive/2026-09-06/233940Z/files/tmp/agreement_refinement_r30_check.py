#!/usr/bin/env python3
"""Exact finite checks for the Wave 30 signature-refinement lemma.

These checks verify the expectation/variance algebra and the deterministic
operator-norm step by exhaustive enumeration in a toy instance.  They also
verify that the A9 wall has coarsest cap 80 but refinement cap 88, illustrating
why the theorem necessarily spends slack.
"""

from __future__ import annotations

import itertools
import numpy as np


def signs(length: int):
    for bits in itertools.product((-1, 1), repeat=length):
        yield np.array(bits, dtype=np.int64)


def toy_variance_check() -> None:
    # Deterministic symmetric zero-diagonal signing.
    A = np.array([
        [0, 1, -1, 1, 1],
        [1, 0, 1, -1, 1],
        [-1, 1, 0, 1, -1],
        [1, -1, 1, 0, 1],
        [1, 1, -1, 1, 0],
    ], dtype=np.int64)
    x0 = np.array([1, -1, 1, 1, -1], dtype=np.int64)
    classes = ((0, 2, 3), (1, 4))
    rs = (2, 2)
    n = len(A)
    K = sum(rs)
    D = np.diag(x0)

    offsets = np.cumsum((0,) + rs)
    outcomes = []
    for labels0 in itertools.product(range(rs[0]), repeat=len(classes[0])):
        for labels1 in itertools.product(range(rs[1]), repeat=len(classes[1])):
            P = np.zeros((n, K), dtype=np.float64)
            for j, (C, labels) in enumerate(zip(classes, (labels0, labels1))):
                for i, label in zip(C, labels):
                    P[i, offsets[j] + label] = 1
            outcomes.append(A @ D @ P)

    empirical_mean = sum(outcomes) / len(outcomes)
    theoretical_mean = np.zeros_like(empirical_mean)
    us = []
    for j, C in enumerate(classes):
        indicator = np.zeros(n)
        indicator[list(C)] = 1
        uj = A @ D @ indicator
        us.append(uj)
        for a in range(rs[j]):
            theoretical_mean[:, offsets[j] + a] = uj / rs[j]
    assert np.allclose(empirical_mean, theoretical_mean)

    empirical_left = np.zeros((n, n))
    empirical_right = np.zeros((K, K))
    for M in outcomes:
        Z = M - theoretical_mean
        empirical_left += Z @ Z.T
        empirical_right += Z.T @ Z
    empirical_left /= len(outcomes)
    empirical_right /= len(outcomes)

    vcols = [x0[i] * A[:, i] for i in range(n)]
    theoretical_left = np.zeros((n, n))
    theoretical_right = np.zeros((K, K))
    for j, C in enumerate(classes):
        rj = rs[j]
        Qj = np.eye(rj) / rj - np.ones((rj, rj)) / (rj * rj)
        theoretical_right[offsets[j]:offsets[j + 1], offsets[j]:offsets[j + 1]] = (
            len(C) * (n - 1) * Qj
        )
        for i in C:
            theoretical_left += (1 - 1 / rj) * np.outer(vcols[i], vcols[i])
    assert np.allclose(empirical_left, theoretical_left)
    assert np.allclose(empirical_right, theoretical_right)

    Csig = max(
        np.dot(sum(e * u for e, u in zip(eps, us)),
               sum(e * u for e, u in zip(eps, us)))
        for eps in itertools.product((-1, 1), repeat=len(classes))
    )
    for M in outcomes:
        Z = M - theoretical_mean
        actual = max(np.dot(M @ z, M @ z) for z in signs(K))
        deterministic_bound = (
            np.sqrt(Csig) + np.sqrt(K) * np.linalg.norm(Z, ord=2)
        ) ** 2
        assert actual <= deterministic_bound + 1e-9

    print("TOY_VARIANCE_IDENTITIES", len(outcomes), "outcomes", "C_sig", int(Csig))


def a9_slack_check() -> None:
    A9 = np.array([
        [0,1,-1,1,1,1,1,-1,-1], [1,0,1,-1,-1,-1,1,-1,-1],
        [-1,1,0,1,1,1,1,1,-1], [1,-1,1,0,1,1,1,-1,1],
        [1,-1,1,1,0,1,-1,1,-1], [1,-1,1,1,1,0,-1,-1,-1],
        [1,1,1,1,-1,-1,0,1,1], [-1,-1,1,-1,1,-1,1,0,-1],
        [-1,-1,-1,1,-1,-1,1,-1,0],
    ], dtype=np.int64)

    def x_of(mask: int) -> np.ndarray:
        return np.array([1] + [-1 if (mask >> (i - 1)) & 1 else 1
                              for i in range(1, 9)], dtype=np.int64)

    family = [x_of(mask) for mask in (1, 99, 99)]
    base = family[0]
    signature_classes: dict[tuple[int, ...], list[int]] = {}
    for i in range(9):
        key = tuple(int(x[i] * base[i]) for x in family[1:])
        signature_classes.setdefault(key, []).append(i)
    assert sorted(map(len, signature_classes.values())) == [3, 6]

    def row(x: np.ndarray) -> int:
        y = A9 @ x
        return int(y @ y)

    classes = tuple(signature_classes.values())
    coarsest = []
    for eps in itertools.product((-1, 1), repeat=len(classes)):
        x = base.copy()
        for e, C in zip(eps, classes):
            x[C] *= e
        coarsest.append(row(x))
    assert max(coarsest) == 80

    partition = ((7,), (0, 3), (1, 5), (2, 6), (4, 8))
    refined = []
    for eps in itertools.product((-1, 1), repeat=len(partition)):
        x = base.copy()
        for e, B in zip(eps, partition):
            x[list(B)] *= e
        refined.append(row(x))
    assert max(refined) == 88
    print("A9_SLACK", "coarsest_cap", max(coarsest), "displayed_refinement_cap", max(refined))


if __name__ == "__main__":
    toy_variance_check()
    a9_slack_check()
