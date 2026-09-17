#!/usr/bin/env python3
"""Finite audit of the Wave 57 incidence-retention lemma and partition wall."""

from fractions import Fraction
from itertools import combinations
from math import comb
import random


def selectors(n, m):
    return [sum(1 << i for i in S) for S in combinations(range(n), m)]


def ov(x, y):
    return bin(x & y).count("1")


def main():
    # This instance has a balanced partition whose expected common-core
    # retention is already below lambda_2.
    n, m, ell = 12, 7, 3
    omega = selectors(n, m)
    N = len(omega)
    B = 36
    assert N % B == 0
    r = N // B

    d_ell = comb(m, ell) * comb(n - ell, m - ell)
    h_ell = Fraction(1, comb(n - ell, m - ell))
    lambda2 = Fraction(
        ell * (ell - 1) * (n - m) * (n - m - 1),
        m * (m - 1) * (n - ell) * (n - ell - 1),
    )
    assert lambda2 == Fraction(
        comb(m - 2, ell - 2) * comb(n - ell - 2, m - ell), d_ell
    )

    # Shell degrees of the overlap-threshold graph (self excluded).
    A = {
        j: comb(m, j) * comb(n - m, m - j)
        for j in range(m + 1)
        if 0 <= m - j <= n - m
    }
    D = {s: sum(A.get(j, 0) for j in range(s, m)) for s in range(ell, m)}
    R = {s: Fraction(D[s] * comb(s, ell), d_ell) for s in D}
    s_star = max(R, key=R.get)

    # Precompute the only pair data used below.
    pair_weight = {}
    pair_overlap = {}
    for i in range(N):
        for j in range(i + 1, N):
            q = ov(omega[i], omega[j])
            pair_overlap[i, j] = q
            pair_weight[i, j] = comb(q, ell) if q >= ell else 0

    def partition_statistics(seed):
        ids = list(range(N))
        random.Random(seed).shuffle(ids)
        blocks = [ids[i : i + r] for i in range(0, N, r)]
        # W and E count ordered distinct pairs.
        W = 0
        E = 0
        shell_counts = {j: 0 for j in range(m)}
        for block in blocks:
            for u, v in combinations(block, 2):
                key = (u, v) if u < v else (v, u)
                q = pair_overlap[key]
                W += 2 * pair_weight[key]
                shell_counts[q] = shell_counts.get(q, 0) + 2
                if q >= s_star:
                    E += 2
        P = h_ell + Fraction(W, N * d_ell)
        Pi = Fraction(E, N * D[s_star])
        histogram_P = h_ell + sum(
            Fraction(comb(j, ell) * count, N * d_ell)
            for j, count in shell_counts.items()
            if j >= ell
        )
        assert P == histogram_P
        # Exact threshold-retention implication P >= h + R_s Pi.
        assert P >= h_ell + R[s_star] * Pi
        return P, Pi

    # Pick a deterministic low-retention partition. Existence itself follows
    # analytically from the random-partition expectation in the memo.
    trials = [(*partition_statistics(seed), seed) for seed in range(80)]
    P, Pi, seed = min(trials, key=lambda x: x[0])

    q_partition = Fraction(r - 1, N - 1)
    expected_P = h_ell + (1 - h_ell) * q_partition
    theta_critical = (lambda2 - h_ell) / R[s_star]

    assert expected_P < lambda2
    assert P < lambda2
    assert Pi < theta_critical

    # If Pi exceeded this exact threshold by eps/R, the lemma would force
    # P-lambda_2 >= eps. Check the algebra for an arbitrary rational eps.
    eps = Fraction(1, 10_000)
    sufficient_pi = (lambda2 - h_ell + eps) / R[s_star]
    assert h_ell + R[s_star] * sufficient_pi - lambda2 == eps

    print(f"n={n}, m={m}, ell={ell}, N={N}, blocks={B}, block_size={r}")
    print(f"lambda2={float(lambda2):.9f}, h={float(h_ell):.9f}")
    print(f"best threshold s={s_star}, R_s={float(R[s_star]):.9f}")
    print(f"random-partition E[P]={float(expected_P):.9f}")
    print(f"selected seed={seed}, P={float(P):.9f}, P-lambda2={float(P-lambda2):.9f}")
    print(f"Pi_s={float(Pi):.9f}, critical Pi_s={float(theta_critical):.9f}")
    print("all exact identity and implication checks passed")


if __name__ == "__main__":
    main()
