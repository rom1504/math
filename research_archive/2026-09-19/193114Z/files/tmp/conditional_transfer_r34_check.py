#!/usr/bin/env python3
"""Exact/finite checks for Wave 34 conditional witness transfer.

The response tables are enumerated with integer arithmetic.  Linear-program
values are used only as an additional audit; the decisive A9 certificates
are explicit integer primal/dual certificates checked below.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product
from math import comb

import numpy as np
from scipy.optimize import linprog


A6 = np.array([
    [0, 1, 1, 1, 1, 1],
    [1, 0, -1, -1, 1, 1],
    [1, -1, 0, 1, -1, 1],
    [1, -1, 1, 0, 1, -1],
    [1, 1, -1, 1, 0, -1],
    [1, 1, 1, -1, -1, 0],
], dtype=np.int64)

A9 = np.array([
    [0, 1, -1, 1, 1, 1, 1, -1, -1],
    [1, 0, 1, -1, -1, -1, 1, -1, -1],
    [-1, 1, 0, 1, 1, 1, 1, 1, -1],
    [1, -1, 1, 0, 1, 1, 1, -1, 1],
    [1, -1, 1, 1, 0, 1, -1, 1, -1],
    [1, -1, 1, 1, 1, 0, -1, -1, -1],
    [1, 1, 1, 1, -1, -1, 0, 1, 1],
    [-1, -1, 1, -1, 1, -1, 1, 0, -1],
    [-1, -1, -1, 1, -1, -1, 1, -1, 0],
], dtype=np.int64)


def states(n: int):
    for sigma in (-1, 1):
        for tail in product((-1, 1), repeat=n - 1):
            yield sigma, np.asarray((1,) + tail, dtype=np.int64)


def qnorm(a: np.ndarray) -> int:
    return max(int(sigma * x @ a @ x) for sigma, x in states(len(a)))


def response_table(a: np.ndarray, u: tuple[int, ...]):
    """Return edge vectors d, outside response Z, labels, energy and deficit."""
    n = len(a)
    outside = tuple(i for i in range(n) if i not in u)
    b = len(u)
    edges = tuple(combinations(range(b), 2))
    auo = a[np.ix_(u, outside)]
    aoo = a[np.ix_(outside, outside)]
    labels, drows, zrows = [], [], []
    for sigma, x in states(b):
        labels.append((sigma, x))
        drows.append([sigma * x[i] * x[j] for i, j in edges])
        zrows.append(max(
            int(sigma * (y @ aoo @ y + 2 * x @ auo @ y))
            for y in map(np.asarray, product((-1, 1), repeat=len(outside)))
        ))
    d = np.asarray(drows, dtype=np.int64)
    z = np.asarray(zrows, dtype=np.int64)
    actual = np.asarray([a[u[i], u[j]] for i, j in edges], dtype=np.int64)
    energy = 2 * d @ actual
    qb = int(energy.max())
    deficit = qb - energy
    return edges, d, z, labels, actual, energy, deficit


def continuous_value(d: np.ndarray, z: np.ndarray, mask: np.ndarray | None = None):
    if mask is None:
        mask = np.ones(len(z), dtype=bool)
    dd, zz = d[mask].astype(float), z[mask].astype(float)
    m = d.shape[1]
    primal = linprog(
        np.r_[np.zeros(m), 1.0],
        A_ub=np.c_[2.0 * dd, -np.ones(len(dd))],
        b_ub=-zz,
        bounds=[(-1.0, 1.0)] * m + [(None, None)],
        method="highs",
    )
    assert primal.success
    return float(primal.fun)


def integral_value(d: np.ndarray, z: np.ndarray, mask: np.ndarray | None = None):
    if mask is None:
        mask = np.ones(len(z), dtype=bool)
    beta = np.asarray(list(product((-1, 1), repeat=d.shape[1])), dtype=np.int16)
    score = z[mask].astype(np.int16)[None, :] + 2 * beta @ d[mask].astype(np.int16).T
    return int(score.max(axis=1).min())


def check_explicit_terminal_wall():
    """Integer certificates for the 3+6 A9 witness-migration wall."""
    u = (3, 4, 5, 6, 7, 8)
    edges, d, z, labels, actual, energy, deficit = response_table(A9, u)
    assert qnorm(A9) == 24
    assert int(energy.max()) == 14
    good = deficit == 0
    assert int(good.sum()) == 4

    # Same local spin, opposite orientations.  Their d-vectors cancel, their
    # external responses average to 24, and both lie outside delta <= 8.
    pair = (28, 60)
    assert np.array_equal(d[pair[0]], -d[pair[1]])
    assert tuple(map(int, deficit[list(pair)])) == (16, 12)
    assert tuple(map(int, z[list(pair)])) == (26, 22)
    assert tuple(map(int, (z + energy)[list(pair)])) == (24, 24)
    assert int(z[list(pair)].sum()) == 48

    # An exact fractional primal exposes no favorable state.  Edge order is
    # lexicographic on the six local coordinates.
    zstar = np.asarray((-1, 1, 1, 1, -1, 1, 0, 0, 1, -1, 1, -1, -1, 0, 0),
                       dtype=np.int64)
    score = z + 2 * d @ zstar
    assert int(score.max()) == 24
    assert int(score[good].max()) == 16

    # Hence V_all = V_bad = 24 in both continuous and integral games, even
    # when favorable means all states through child deficit eight.
    for cutoff in (0, 4, 8):
        bad = deficit > cutoff
        assert continuous_value(d, z, bad) == 24.0
        assert integral_value(d, z, bad) == 24
    assert continuous_value(d, z) == 24.0
    assert integral_value(d, z) == 24

    # Exact indicator-penalty formula V(tau)=24-tau for 0 <= tau <= 8.
    # The bad pair is the lower certificate.  zstar is the upper certificate:
    # all bad scores lose tau and all good scores are at most 16.
    for tau in map(Fraction, range(9)):
        scores = [Fraction(int(z[i] + 2 * d[i] @ zstar))
                  - tau * int(not good[i]) for i in range(len(z))]
        assert max(scores) == 24 - tau

    # On every switched/complemented A6 replacement, the global maximum is
    # attained only outside the original B-ground fiber, with gap >= 8.
    cedges = np.asarray([A6[i, j] for i, j in edges], dtype=np.int64)
    norms, gaps = [], []
    for h in range(len(labels)):
        beta = cedges * d[h]
        scores = z + 2 * d @ beta
        norms.append(int(scores.max()))
        gaps.append(int(scores.max() - scores[good].max()))
        assert int(scores[~good].max()) == int(scores.max())
        assert not np.any((scores == scores.max()) & good)
    assert {v: norms.count(v) for v in set(norms)} == {28: 6, 32: 20, 36: 26, 40: 12}
    assert min(gaps) == 8

    return {"V": 24, "good_count": 4, "good_exposure": 16,
            "replacement_norm_hist": {v: norms.count(v) for v in sorted(set(norms))}}


def check_same_orientation_shield():
    """A9 bad-fiber certificate supported in one absolute orientation."""
    u = (0, 1, 2, 6)
    edges, d, z, labels, actual, energy, deficit = response_table(A9, u)
    pair = (11, 14)
    assert [labels[i][0] for i in pair] == [1, 1]
    assert tuple(map(int, deficit[list(pair)])) == (12, 12)
    assert tuple(map(int, energy[list(pair)])) == (-4, -4)
    assert tuple(map(int, z[list(pair)])) == (28, 28)
    mean_d = d[list(pair)].sum(axis=0) // 2
    assert tuple(map(int, mean_d)) == (0, 1, 0, 0, -1, 0)
    # Average replacement score is 28+2(z_02-z_16) >= 24 on the cube.
    assert 28 - 2 * int(np.abs(mean_d).sum()) == 24
    for cutoff in (0, 4, 8):
        bad = deficit > cutoff
        assert continuous_value(d, z, bad) == 24.0
        assert integral_value(d, z, bad) == 24
    return {"orientation": 1, "deficits": (12, 12), "J": 24}


def check_a6_uniform_finite_wall():
    """Every A6 proper selector of sizes 3--5 has a bad-only certificate."""
    assert qnorm(A6) == 10
    counts = {}
    for b in (3, 4, 5):
        count = 0
        for u in combinations(range(6), b):
            edges, d, z, labels, actual, energy, deficit = response_table(A6, u)
            half = 2 ** (b - 1)
            found = False
            for i in range(half):
                j = i + half
                if deficit[i] > 0 and deficit[j] > 0 and z[i] + z[j] == 20:
                    assert np.array_equal(d[i], -d[j])
                    found = True
                    break
            assert found
            bad = deficit > 0
            assert continuous_value(d, z, bad) == 10.0
            assert integral_value(d, z, bad) == 10
            count += 1
        counts[b] = count
    assert counts == {3: 20, 4: 15, 5: 6}
    return counts


def main():
    terminal = check_explicit_terminal_wall()
    shield = check_same_orientation_shield()
    a6 = check_a6_uniform_finite_wall()
    print("A9 terminal", terminal)
    print("A9 one-orientation shield", shield)
    print("A6 all-selector finite wall", a6)
    print("conditional_transfer_r34_check: PASS")


if __name__ == "__main__":
    main()
