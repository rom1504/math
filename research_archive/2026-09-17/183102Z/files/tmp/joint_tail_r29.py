#!/usr/bin/env python3
"""Exact finite census for joint H_S alignment and row square."""

from fractions import Fraction
from itertools import combinations, product
import sys

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from envelope_block_cover_r27 import A6, A8, A9  # noqa: E402


def spins(n):
    return np.array([(1,) + z for z in product((-1, 1), repeat=n - 1)], dtype=np.int64)


def qnorm(a):
    x = spins(len(a))
    return int(np.max(np.abs(np.einsum("bi,ij,bj->b", x, a, x))))


def audit(a, name):
    n = len(a)
    xx = spins(n)
    r = np.einsum("bi,ij,bj->b", xx, a @ a, xx).astype(int)
    q = qnorm(a)
    print(name, "n", n, "q", q, "row values", sorted(set(map(int, r))))
    candidates = []
    for m in range(2, n):
        p2 = Fraction(m * (m - 1), n * (n - 1))
        # Work with H multiplied by den=p2.denominator for exact integer values.
        den = p2.denominator
        num = p2.numerator
        for S in combinations(range(n), m):
            mask = np.zeros(n, dtype=np.int64)
            mask[list(S)] = 1
            hnum_mat = a * (den * np.outer(mask, mask) - num)
            hnum = np.abs(np.einsum("bi,ij,bj->b", xx, hnum_mat, xx)).astype(int)
            child = a[np.ix_(S, S)]
            qs = qnorm(child)
            # Y uses irrational p^(3/2); retain float only for selecting actual threshold.
            y = qs - (m / n) ** 1.5 * q
            if y <= 1e-12:
                continue
            hit = hnum / den + 1e-12 >= y
            minrow = int(r[hit].min()) if np.any(hit) else None
            hitprob = Fraction(int(hit.sum()), len(xx))
            # Strict separation: largest row cap below every actual-threshold hit.
            candidates.append((minrow, hitprob, y, m, S, sorted(set(map(int, hnum[hit]))), den))
    candidates.sort(reverse=True)
    for z in candidates[:8]:
        minrow, prob, y, m, S, hv, den = z
        print(" candidate", "m", m, "S", S, "Y", round(y, 9), "P", prob,
              "minRhit", minrow, "hit |H| nums/den", hv, "/", den)


def main():
    for name, a in (("A6", A6), ("A8", A8), ("A9", A9)):
        audit(np.asarray(a, dtype=np.int64), name)

    # Exact A8 joint-tail wall used in the memo.
    a = np.asarray(A8, dtype=np.int64)
    n, m, ss = 8, 4, (2, 3, 5, 7)
    xx = spins(n)
    mask = np.zeros(n, dtype=np.int64)
    mask[list(ss)] = 1
    # Since p2=3/14, 7*x'Hx is integral (the quadratic double-counts edges).
    h14 = np.abs(np.einsum(
        "bi,ij,bj->b", xx, a * (14 * np.outer(mask, mask) - 3), xx
    )).astype(int)
    assert np.all(h14 % 2 == 0)
    h7 = h14 // 2
    row = np.einsum("bi,ij,bj->b", xx, a @ a, xx).astype(int)
    assert qnorm(a) == 20 and qnorm(a[np.ix_(ss, ss)]) == 12
    # 34/7 < 12-5sqrt(2) < 40/7, checked by exact squared inequalities.
    assert 2500 > 2450 and 1936 < 2450
    hit = h7 >= 40
    assert int(hit.sum()) == 26
    assert {int(k): int(v) for k, v in zip(*np.unique(row[hit], return_counts=True))} == {
        56: 8, 64: 12, 72: 6
    }
    assert not np.any(hit & (row <= 40))
    assert int(np.sum(row <= 40)) == 30
    assert Fraction(int(row[hit].sum()), int(hit.sum())) == Fraction(824, 13)
    assert Fraction(int(row.sum()), len(row)) == 56
    print(
        "PASS exact A8 wall: P(hit)=13/64, P(R<=40)=15/64, joint=0,",
        "E[R|hit]=824/13 > E[R]=56",
    )


if __name__ == "__main__":
    main()
