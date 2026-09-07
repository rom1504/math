#!/usr/bin/env python3
"""Independent exact audit of the q7 -> q6 -> q5 causal suffix wall."""

from itertools import product

import numpy as np


A7 = np.array([
    [0, 1, 1, 1, 1, 1, 1],
    [1, 0, 1, 1, -1, -1, 1],
    [1, 1, 0, 1, -1, 1, -1],
    [1, 1, 1, 0, 1, -1, -1],
    [1, -1, -1, 1, 0, -1, -1],
    [1, -1, 1, -1, -1, 0, -1],
    [1, 1, -1, -1, -1, -1, 0],
], dtype=np.int64)


def q(a):
    best = 0
    for tail in product((-1, 1), repeat=len(a) - 1):
        x = np.array((1,) + tail, dtype=np.int64)
        best = max(best, abs(int(x @ a @ x)))
    return best


def principal(a, removed):
    keep = [i for i in range(len(a)) if i not in set(removed)]
    return a[np.ix_(keep, keep)]


def main():
    a6 = principal(A7, (0,))
    assert q(A7) == 18
    assert q(a6) == 10
    child_norms = tuple(q(principal(A7, (0, j))) for j in range(1, 7))
    assert child_norms == (8,) * 6

    # Both selected edges are deterministic instances of field-proportional
    # deletion after orienting an absolute ground positively.
    x7 = np.array((1, -1, -1, -1, -1, -1, -1), dtype=np.int64)
    assert int(x7 @ A7 @ x7) == -18
    rows7 = -x7 * (A7 @ x7)
    assert tuple(map(int, rows7)) == (6, 0, 0, 0, 4, 4, 4)
    assert rows7[0] == len(A7) - 1

    x6 = np.array((1, -1, -1, -1, -1, -1), dtype=np.int64)
    assert int(x6 @ a6 @ x6) == -10
    rows6 = -x6 * (a6 @ x6)
    assert tuple(map(int, rows6)) == (1, 1, 1, 1, 1, 5)
    assert rows6[-1] == len(a6) - 1

    # x1 = 10 - 108 sqrt(42)/49 is strictly negative.
    assert 490**2 < 108**2 * 42

    # x2 > 0 is equivalent to 9(6 sqrt(6)-5 sqrt(5)) > 7 sqrt(7).
    # Squaring once gives 27278 > 4860 sqrt(30), and both sides are positive.
    assert 27278**2 > 4860**2 * 30

    # x1+x2 = 8 - 90 sqrt(35)/49 is strictly negative.
    assert 392**2 < 90**2 * 35

    c = 18 / (7 ** 1.5)
    x1 = c * (7 ** 1.5 - 6 ** 1.5) - 8
    x2 = c * (6 ** 1.5 - 5 ** 1.5) - 2
    print("norm chain:", q(A7), q(a6), child_norms)
    print("oriented row profiles:", tuple(rows7), tuple(rows6))
    print("centered atoms:", x1, x2, "sum:", x1 + x2)
    assert x1 < 0 < x2 and x1 + x2 < 0
    print("root cut passes, proper suffix cut fails: PASS")


if __name__ == "__main__":
    main()
