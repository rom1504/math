#!/usr/bin/env python3
"""Independent exact audit of the Wave 11 order-six replacement wall."""

from itertools import product


A = (
    (0, 1, 1, 1, 1, 1),
    (1, 0, -1, -1, 1, 1),
    (1, -1, 0, 1, -1, 1),
    (1, -1, 1, 0, 1, -1),
    (1, 1, -1, 1, 0, -1),
    (1, 1, 1, -1, -1, 0),
)

DPRIME = (
    (0, -1, -1, 1, 1),
    (-1, 0, 1, 1, -1),
    (-1, 1, 0, -1, 1),
    (1, 1, -1, 0, -1),
    (1, -1, 1, -1, 0),
)


def energy(a, x):
    return sum(a[i][j] * x[i] * x[j] for i in range(len(a)) for j in range(len(a)))


def extrema(a):
    vals = [(energy(a, (1,) + z), (1,) + z) for z in product((-1, 1), repeat=len(a) - 1)]
    p = max(v for v, _ in vals)
    n = -min(v for v, _ in vals)
    return p, n, [x for v, x in vals if v == p], [x for v, x in vals if v == -n]


def q_of(a):
    p, n, _, _ = extrema(a)
    return max(p, n)


def global_q(n):
    edges = [(i, j) for i in range(n) for j in range(i + 1, n)]
    best = n * (n - 1)
    for bits in range(1 << len(edges)):
        a = [[0] * n for _ in range(n)]
        for k, (i, j) in enumerate(edges):
            z = 1 if (bits >> k) & 1 else -1
            a[i][j] = a[j][i] = z
        best = min(best, q_of(tuple(tuple(row) for row in a)))
    return best


def layer(c, d):
    q = q_of(d)
    best = 0
    for y in product((-1, 1), repeat=len(d)):
        h = energy(d, y)
        exposure = 2 * abs(sum(ci * yi for ci, yi in zip(c, y)))
        best = max(best, exposure - (q - h), exposure - (q + h))
    return max(0, best)


def hybrid(c, d):
    n = len(d) + 1
    out = [[0] * n for _ in range(n)]
    for i in range(len(d)):
        for j in range(len(d)):
            out[i][j] = d[i][j]
        out[i][-1] = out[-1][i] = c[i]
    return tuple(tuple(row) for row in out)


if __name__ == "__main__":
    d = tuple(tuple(A[i][j] for j in range(5)) for i in range(5))
    c = tuple(A[5][i] for i in range(5))
    p, n, ps, ns = extrema(A)
    assert (p, n, q_of(A)) == (10, 10, 10)
    assert global_q(5) == 8
    assert global_q(6) == 10
    assert q_of(d) == q_of(DPRIME) == 8
    assert c == (1, 1, 1, -1, -1)
    assert layer(c, d) == 2
    assert layer(c, DPRIME) == 10
    assert q_of(hybrid(c, DPRIME)) == 18
    target_p = (1, 1, 1, -1, -1, 1)
    target_n = (-1, -1, -1, 1, 1, 1)
    assert target_p in ps or tuple(-z for z in target_p) in ps
    assert target_n in ns or tuple(-z for z in target_n) in ns
    print("q5=8, q6=10; endpoint split and layer wall 2 -> 10; hybrid Q=18: PASS")
