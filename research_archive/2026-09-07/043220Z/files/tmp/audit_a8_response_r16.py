#!/usr/bin/env python3
"""Independent exact audit of the A8 response/endpoint wall."""

from itertools import product

import numpy as np

from check_bicriteria_recurrence_r15 import all_minimizers, qnorm


A8 = np.array(
    [
        [0, 1, 1, 1, 1, 1, 1, 1],
        [1, 0, 1, -1, 1, 1, -1, -1],
        [1, 1, 0, 1, -1, 1, -1, -1],
        [1, -1, 1, 0, -1, -1, -1, 1],
        [1, 1, -1, -1, 0, -1, 1, -1],
        [1, 1, 1, -1, -1, 0, 1, 1],
        [1, -1, -1, -1, 1, 1, 0, 1],
        [1, -1, -1, 1, -1, 1, 1, 0],
    ],
    dtype=np.int64,
)

BLOCKS = ((0, 3, 4, 5), (1, 2, 6, 7))


def energies(M):
    n = len(M)
    spins = np.array([(1,) + t for t in product((-1, 1), repeat=n - 1)], dtype=np.int64)
    return spins, np.einsum("bi,ij,bj->b", spins, M, spins)


def main():
    assert qnorm(A8) == 20
    C = A8.copy()
    block_norms = []
    for block in BLOCKS:
        D = A8[np.ix_(block, block)]
        block_norms.append(qnorm(D))
        C[np.ix_(block, block)] = 0
    assert block_norms == [12, 12]
    assert qnorm(C) == 16

    q4, minimizers = all_minimizers(4)
    assert q4 == 8
    best = 10**9
    count = 0
    for G1 in minimizers:
        for G2 in minimizers:
            H = C.copy()
            H[np.ix_(BLOCKS[0], BLOCKS[0])] = G1
            H[np.ix_(BLOCKS[1], BLOCKS[1])] = G2
            value = qnorm(H)
            if value < best:
                best, count = value, 1
            elif value == best:
                count += 1
    assert best == 20

    S = sum(block_norms)
    B = 2 * q4
    X = S - B
    phi_original = qnorm(A8) - S
    phi_best = best - B
    assert (S, B, X, phi_original, phi_best, phi_best - phi_original) == (24, 16, 8, -4, 4, 8)
    print("Q(A), Q(C), block norms, B, X:", qnorm(A8), qnorm(C), block_norms, B, X)
    print("best pure hybrid, number attaining:", best, count)
    print("common-mosaic responses original/best/increase:", phi_original, phi_best, phi_best - phi_original)


if __name__ == "__main__":
    main()
