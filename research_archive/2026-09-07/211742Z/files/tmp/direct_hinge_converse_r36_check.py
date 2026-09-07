#!/usr/bin/env python3
"""Exact algebra checks for the Wave 36 rare-coverage converse."""

from __future__ import annotations

import itertools
import sys

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from envelope_block_cover_r27 import A6, A8, A9


def spins(n: int):
    for tail in itertools.product((-1, 1), repeat=n - 1):
        yield np.array((1,) + tail, dtype=np.int64)


def audit(A: np.ndarray, name: str, m: int) -> None:
    n = len(A)
    # n^2(n-1) is a common denominator for the exact decomposition.
    common_den = n * n * (n - 1)
    q = max(abs(int(x @ A @ x)) for x in spins(n))
    selectors = list(itertools.combinations(range(n), m))

    checked = 0
    for x in spins(n):
        for sigma in (-1, 1):
            B = sigma * (x[:, None] * A * x[None, :])
            r = B @ np.ones(n, dtype=np.int64)
            row = int(r @ r)
            assert row == int((A @ x) @ (A @ x))
            assert int(np.sum(B * B)) == n * (n - 1)
            assert abs(int(np.ones(n, dtype=np.int64) @ B @ np.ones(n, dtype=np.int64))) <= q
            assert abs(np.linalg.norm(B, 2) - np.linalg.norm(A, 2)) < 1e-9

            total = int(B.sum())
            for bits in itertools.product((0, 1), repeat=n):
                eta = np.array(bits, dtype=np.int64)
                znum = n * eta - m
                lhs_num = (
                    int(eta @ B @ eta) * common_den
                    - m * (m - 1) * n * total
                )
                rhs_num = (
                    m * (n - m) * total
                    + 2 * m * (n - 1) * int(r @ znum)
                    + (n - 1) * int(znum @ B @ znum)
                )
                assert lhs_num == rhs_num

            for S in selectors:
                eta = np.zeros(n, dtype=np.int64)
                eta[list(S)] = 1
                xloss_num = (
                    int(eta @ B @ eta) * n * (n - 1)
                    - m * (m - 1) * total
                )
                direct_num = (
                    sum(int(B[i, j]) for i in S for j in S) * n * (n - 1)
                    - m * (m - 1) * total
                )
                assert xloss_num == direct_num
            checked += 1

    print(name, "n", n, "m", m, "cuts", checked, "selectors", len(selectors), "q", q)


def main() -> None:
    audit(A6, "A6", 3)
    audit(A8, "A8", 5)
    audit(A9, "A9", 6)
    print("PASS: exact C36.1--C36.2 and slice-payoff algebra")


if __name__ == "__main__":
    main()
