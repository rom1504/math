#!/usr/bin/env python3
"""Exact audit of selector-conditioned outside-block switching."""

from collections import defaultdict
from fractions import Fraction
from itertools import combinations, product
from math import comb, log
import random

import numpy as np


def quadratic(a, x):
    return int(x @ a @ x)


def projective_spins(n):
    for tail in product((-1, 1), repeat=n - 1):
        yield np.array((1,) + tail, dtype=int)


def qnorm(a):
    return max(abs(quadratic(a, x)) for x in projective_spins(len(a)))


def flip(x, u):
    ans = x.copy()
    ans[list(u)] *= -1
    return ans


def entropy(probabilities):
    return -sum(p * log(p) for p in probabilities if p)


def trial(n, m, k, rng):
    a = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i + 1, n):
            a[i, j] = a[j, i] = rng.choice((-1, 1))

    q = qnorm(a)
    best = None
    for x in projective_spins(n):
        raw = quadratic(a, x)
        for tau in (-1, 1):
            if tau * raw == q:
                best = (tau, x)
                break
        if best is not None:
            break
    tau, x = best
    p2 = Fraction(m * (m - 1), n * (n - 1))

    selectors = list(combinations(range(n), m))
    # Any nontrivial deterministic family suffices for the averaging audit.
    family = [s for idx, s in enumerate(selectors) if idx % 3 != 1]
    tsize = n - m

    cmat = np.diag(x) @ (a @ a) @ np.diag(x)
    r0 = int(x @ (a @ a) @ x)
    energy0 = tau * quadratic(a, x)
    assert energy0 == q

    u_marginal = defaultdict(Fraction)
    row_joint = Fraction(0)
    deficit_joint = Fraction(0)
    conditional_rows = []

    for s_tuple in family:
        sset = set(s_tuple)
        t_tuple = tuple(i for i in range(n) if i not in sset)
        us = list(combinations(t_tuple, k))

        # Exact block formulas for C=diag(x) A^2 diag(x).
        avec = np.zeros(n, dtype=int)
        avec[list(s_tuple)] = 1
        tvec = 1 - avec
        css = int(avec @ cmat @ avec)
        cst = int(avec @ cmat @ tvec)
        ctt = int(tvec @ cmat @ tvec)
        alpha = Fraction(tsize - 2 * k, tsize)
        beta = Fraction(
            tsize * (tsize - 1) - 4 * k * (tsize - k),
            tsize * (tsize - 1),
        )
        predicted_row = (
            css + 2 * alpha * cst + beta * ctt
            + (1 - beta) * tsize * (n - 1)
        )

        rows = []
        for u in us:
            xu = flip(x, u)
            rowu = int(xu @ (a @ a) @ xu)
            deltau = q - tau * quadratic(a, xu)
            assert deltau >= 0
            rows.append(rowu)
            row_joint += Fraction(rowu, len(family) * len(us))
            deficit_joint += Fraction(deltau, len(family) * len(us))
            u_marginal[u] += Fraction(1, len(family) * len(us))

            # For U disjoint from S the local energy and local deficit agree.
            loc0 = tau * quadratic(a[np.ix_(s_tuple, s_tuple)], x[list(s_tuple)])
            locu = tau * quadratic(a[np.ix_(s_tuple, s_tuple)], xu[list(s_tuple)])
            assert locu == loc0
            # Delta(d)=0 here, so the bare-loss change is -p2 Delta(d^U).
            assert -p2 * deltau <= 0

        actual_row = Fraction(sum(rows), len(rows))
        assert actual_row == predicted_row, (actual_row, predicted_row)
        conditional_rows.append(actual_row)

    assert row_joint == sum(conditional_rows, Fraction(0)) / len(family)

    # D=d^U identifies U because k<n/2.  Compute I(S;U) exactly.
    h_u = entropy(u_marginal.values())
    h_u_given_s = log(comb(tsize, k))
    mutual_information = h_u - h_u_given_s
    information_bound = log(comb(n, k) / comb(tsize, k))
    assert mutual_information <= information_bound + 1e-12

    return float(row_joint), float(deficit_joint), mutual_information, information_bound


def main():
    rng = random.Random(590059)
    records = []
    for n in range(6, 10):
        for m in range(n // 2 + 1, n - 1):
            t = n - m
            for k in range(1, min(t, (n - 1) // 2) + 1):
                for _ in range(5):
                    records.append(trial(n, m, k, rng))
    print("outside-block row, deficit, fibre preservation, and information checks passed")
    print("audited instances:", len(records))


if __name__ == "__main__":
    main()
