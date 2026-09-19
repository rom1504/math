#!/usr/bin/env python3
"""Finite numerical audits of the Wave 21 fixed-slice mgf inequality."""

from itertools import combinations, product
from math import comb, exp, log
import random


def logmeanexp(values):
    top = max(values)
    return top + log(sum(exp(v - top) for v in values) / len(values))


def op_norm_symmetric(a):
    # Power iteration is not a certified upper bound, so use the Frobenius
    # norm in the theorem audit; ||A||_op <= ||A||_F.
    return sum(a[i][j] ** 2 for i in range(len(a)) for j in range(len(a))) ** 0.5


def one_audit(a, x, sigma, m, lam):
    n = len(a)
    w = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            w[i][j] = sigma * a[i][j] * x[i] * x[j]
    rows = [sum(w[i][j] for j in range(n) if j != i) for i in range(n)]
    r2 = sum(z * z for z in rows)
    rinf = max(map(abs, rows))
    W = sum(w[i][j] for i in range(n) for j in range(i + 1, n))
    subsets = list(combinations(range(n), m))
    vals = [2 * sum(w[i][j] for i, j in combinations(s, 2)) for s in subsets]
    p = m / n
    p2 = m * (m - 1) / (n * (n - 1))
    centered = [lam * (v - 2 * p2 * W) for v in vals]
    lhs = logmeanexp(centered)

    prob = comb(n, m) * p**m * (1 - p) ** (n - m)
    chi = -log(prob)
    eps = p * p - p2
    # For a signing, the maximum rectangular squared Frobenius norm is at
    # most floor(n^2/4).  Frobenius is a rigorous (looser) op-norm ceiling.
    fcross = n * n // 4
    op = op_norm_symmetric(w)
    assert lam < 3 / (4 * p * rinf) if rinf else True
    assert lam < 1 / (2 * op) if op else True
    rhs = (
        chi
        + lam * eps * abs(2 * W)
        + 4 * lam * lam * p * p * p * (1 - p) * r2
        / (1 - 4 * lam * p * rinf / 3)
        + lam * lam * fcross / (1 - 4 * lam * lam * op * op)
    )
    assert lhs <= rhs + 1e-11, (n, m, lam, lhs, rhs)
    return rhs - lhs


def main():
    rng = random.Random(2101)
    minimum_slack = float("inf")
    tests = 0
    for n in range(4, 9):
        for _ in range(12):
            a = [[0] * n for _ in range(n)]
            for i, j in combinations(range(n), 2):
                a[i][j] = a[j][i] = rng.choice((-1, 1))
            x = [rng.choice((-1, 1)) for _ in range(n)]
            sigma = rng.choice((-1, 1))
            for m in range(2, n):
                # This conservative choice is inside both domains even with
                # the Frobenius op-norm ceiling used by the checker.
                lam = 0.04 / n
                minimum_slack = min(minimum_slack, one_audit(a, x, sigma, m, lam))
                tests += 1
    print(f"fixed-slice mgf audits passed: {tests}; minimum slack={minimum_slack:.12g}")


if __name__ == "__main__":
    main()
