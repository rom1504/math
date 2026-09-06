#!/usr/bin/env python3
"""Exact finite checks for the Wave 40 favorable distance transform."""

from __future__ import annotations

import itertools
import math

import numpy as np


def pspins(n):
    for tail in itertools.product((-1, 1), repeat=n - 1):
        yield np.array((1,) + tail, dtype=np.int64)


def prdist(x, y):
    d = int(np.count_nonzero(x != y))
    return min(d, len(x) - d)


def qnorm(a):
    return max(abs(int(x @ a @ x)) for x in pspins(len(a)))


def grounds(a, selector):
    b = a[np.ix_(selector, selector)]
    ys = list(pspins(len(selector)))
    vals = [abs(int(y @ b @ y)) for y in ys]
    q = max(vals)
    return q, [y for y, val in zip(ys, vals) if val == q], sorted(set(q - val for val in vals))


def cylinder(a, selector, label):
    out = []
    for x in pspins(len(a)):
        r = x[list(selector)]
        if np.array_equal(r, label) or np.array_equal(r, -label):
            out.append(x)
    return out


def cylinder_distance(a, s, y, t, h):
    return min(prdist(x, z) for x in cylinder(a, s, y) for z in cylinder(a, t, h))


def signing_from_bits(n, bits):
    a = np.zeros((n, n), dtype=np.int64)
    for (i, j), value in zip(itertools.combinations(range(n), 2), bits):
        a[i, j] = a[j, i] = value
    return a


def exact_minimum_q(n):
    return min(qnorm(signing_from_bits(n, bits)) for bits in itertools.product((-1, 1), repeat=n * (n - 1) // 2))


def main():
    # This is an exact order-four minimizer: exhaustive enumeration of all 64
    # complete signings gives q_4=8.
    a = np.array(
        [
            [0, -1, -1, -1],
            [-1, 0, -1, -1],
            [-1, -1, 0, 1],
            [-1, -1, 1, 0],
        ],
        dtype=np.int64,
    )
    assert exact_minimum_q(4) == qnorm(a) == 8

    s = (0, 1, 2)
    t = (0, 2, 3)
    qs, fs, deficits_s = grounds(a, s)
    qt, ft, deficits_t = grounds(a, t)
    assert qs == qt == 6
    assert len(fs) == len(ft) == 1
    assert fs[0].tolist() == [1, 1, 1]
    assert ft[0].tolist() == [1, -1, -1]
    assert deficits_s == deficits_t == [0, 4]

    # At t=0, B_(4,3)=3 sqrt(3)-4 lies strictly between 0 and the next
    # deficit 4, so the actual favorable projective fibers are exactly fs,ft.
    b43 = 3 * math.sqrt(3) - 4
    assert 0 < b43 < 4
    assert cylinder_distance(a, s, fs[0], t, ft[0]) == 1

    # Exhaust the distance-transform identity for the two terminal cylinders.
    ys = list(pspins(3))
    def dset(y, fset):
        return min(prdist(y, f) for f in fset)

    hard = min(cylinder_distance(a, s, y, t, h) for y in fs for h in ft)
    soft = min(
        cylinder_distance(a, s, y, t, h) + dset(y, fs) + dset(h, ft)
        for y in ys
        for h in ys
    )
    assert hard == soft == 1

    # Common low-row root scaffold: z restricts favorably on s and is one
    # projective flip from ft on t.  Its soft cost equals the hard optimum.
    z = np.ones(4, dtype=np.int64)
    assert int((a @ z) @ (a @ z)) == 20 <= 2 * 4 * 3
    star_soft = dset(z[list(s)], fs) + dset(z[list(t)], ft)
    assert star_soft == hard == 1

    # With singleton normalized ground measures, the local log normalizers
    # are zero but the global compatibility partition retains e^{-lambda}.
    lam = 2.0
    partition = math.exp(-lam * hard)
    assert abs(-math.log(partition) / lam - 1.0) < 1e-12

    print("PASS q4=8 exact minimizer and t=0 incompatible favorable fibers")
    print("PASS hard distance = soft distance-transform = star scaffold = 1")
    print("PASS singleton log partition retains boundary compatibility cost 1")


if __name__ == "__main__":
    main()
