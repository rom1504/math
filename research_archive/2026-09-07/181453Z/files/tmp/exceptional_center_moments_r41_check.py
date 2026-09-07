#!/usr/bin/env python3
"""Exact checks for the Wave 41 fixed-slice moment audit."""

from __future__ import annotations

import itertools
from collections import Counter
from fractions import Fraction

import numpy as np

from envelope_block_cover_r27 import A6, A8, A9


def spins(n):
    for tail in itertools.product((-1, 1), repeat=n - 1):
        yield np.array((1,) + tail, dtype=np.int64)


def falling(a, j):
    ans = 1
    for k in range(j):
        ans *= a - k
    return ans


def pj(m, n, j):
    if j > m:
        return Fraction(0)
    return Fraction(falling(m, j), falling(n, j))


def qnorm(a):
    return max(abs(int(x @ a @ x)) for x in spins(len(a)))


def empirical(a, z, m):
    vals = []
    for s in itertools.combinations(range(len(a)), m):
        b = a[np.ix_(s, s)]
        c = int(z[list(s)] @ b @ z[list(s)])
        vals.append(c)
    return vals


def formula(a, z, m):
    n = len(a)
    e = int(z @ a @ z)
    row = int((a @ z) @ (a @ z))
    p2, p3, p4 = (pj(m, n, j) for j in (2, 3, 4))
    mean = p2 * e
    second = (
        p4 * e * e
        + 4 * (p3 - p4) * row
        + 2 * n * (n - 1) * (p2 - 2 * p3 + p4)
    )
    variance = second - mean * mean
    return mean, second, variance


def check_all_moments():
    # Exhaust all projective z and all proper slice sizes in the three finite
    # minimizers.  Every comparison is rational/integer exact.
    for a in (A6, A8, A9):
        n = len(a)
        for z in spins(n):
            for m in range(2, n):
                vals = empirical(a, z, m)
                den = len(vals)
                mean = Fraction(sum(vals), den)
                second = Fraction(sum(c * c for c in vals), den)
                fm, fs, fv = formula(a, z, m)
                assert mean == fm
                assert second == fs
                assert fv == second - mean * mean


def a8_wall():
    # Two low-row cuts of the same exact minimizer have the same (E,R), hence
    # the same first two signed-energy moments, but different absolute-deficit
    # tails.  This isolates the uncontrolled Q(A[S])-|C_S| coupling.
    xs = list(spins(8))
    z0, z12 = xs[0], xs[12]
    m = 4
    q = qnorm(A8)
    assert q == 20
    assert (int(z0 @ A8 @ z0), int((A8 @ z0) @ (A8 @ z0))) == (-16, 72)
    assert (int(z12 @ A8 @ z12), int((A8 @ z12) @ (A8 @ z12))) == (-16, 72)
    assert formula(A8, z0, m) == formula(A8, z12, m)
    assert formula(A8, z0, m)[:2] == (Fraction(-24, 7), Fraction(208, 7))

    ys = list(spins(m))
    histograms = []
    signed_third = []
    for z in (z0, z12):
        deficits = []
        cs = []
        for s in itertools.combinations(range(8), m):
            b = A8[np.ix_(s, s)]
            qs = max(abs(int(y @ b @ y)) for y in ys)
            c = int(z[list(s)] @ b @ z[list(s)])
            deficits.append(qs - abs(c))
            cs.append(c)
        histograms.append(Counter(deficits))
        signed_third.append(sum(c ** 3 for c in cs))
    assert histograms[0] == Counter({4: 34, 0: 19, 8: 15, 12: 2})
    assert histograms[1] == Counter({4: 28, 0: 22, 8: 16, 12: 4})
    assert signed_third == [-14592, -13824]

    # B_(8,4)=20(2^(-3/2)-3/14) lies strictly between 0 and the next deficit 4.
    # Squaring avoids a floating comparison: B>0 and B<4 are immediate from
    # 5/sqrt(2)-30/7 in exact form.
    bfloat = 20 * ((m / 8) ** 1.5 - Fraction(m * (m - 1), 8 * 7))
    assert 0 < bfloat < 4
    assert sum(v <= bfloat for v in histograms[0].elements()) == 19
    assert sum(v <= bfloat for v in histograms[1].elements()) == 22
    assert 72 <= 2 * 8 * 7
    return histograms


if __name__ == "__main__":
    check_all_moments()
    h0, h1 = a8_wall()
    print("PASS exact fixed-slice first/second moments on A6,A8,A9")
    print("PASS A8 same-(E,R) absolute-deficit wall:", dict(h0), dict(h1))
