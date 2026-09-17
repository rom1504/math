#!/usr/bin/env python3
"""Exact verifier for Wave 7 two-channel harvesting obstructions.

The default checks the displayed witnesses using integer enumeration.  Pass
--exhaustive-small to enumerate every labelled signing through order six and
all positive/negative endpoint ties.
"""

from __future__ import annotations

import argparse
from functools import lru_cache
from itertools import product


def spins(n):
    """Boolean spins modulo global sign."""
    if n == 0:
        return ((),)
    return tuple(x + (1,) for x in product((-1, 1), repeat=n - 1))


def energy(a, x):
    return sum(a[i][j] * x[i] * x[j] for i in range(len(a)) for j in range(len(a)))


@lru_cache(None)
def extrema(a):
    rows = tuple((energy(a, x), x) for x in spins(len(a)))
    hi = max((v for v, _ in rows), default=0)
    lo = min((v for v, _ in rows), default=0)
    return hi, -lo, tuple(x for v, x in rows if v == hi), tuple(x for v, x in rows if v == lo)


def induced(a, ids):
    return tuple(tuple(a[i][j] for j in ids) for i in ids)


def endpoint_split(a, p, n):
    sides = []
    ratio = tuple(x * y for x, y in zip(p, n))
    for sign in (-1, 1):
        ids = tuple(i for i, z in enumerate(ratio) if z == sign)
        if ids:
            sides.append(ids)
    assert len(sides) == 2
    out = []
    for ids in sides:
        other = tuple(i for i in range(len(a)) if i not in ids)
        child = induced(a, ids)
        pc, nc, _, _ = extrema(child)
        q = max(pc, nc)
        state = tuple(p[i] for i in ids)
        h = energy(child, state)
        fields = tuple(sum(a[j][i] * p[i] for i in ids) for j in other)
        ell = sum(abs(z) for z in fields)
        mu = (2 * ell - (q - h), 2 * ell - (q + h))
        out.append((ids, child, pc, nc, h, ell, mu, abs(pc - nc)))
    return tuple(out)


@lru_cache(None)
def liberal_path(a):
    """Best one-path sum, optimizing endpoint ties and one orientation/edge."""
    if len(a) <= 1:
        return 0, ()
    _, _, ps, ns = extrema(a)
    best = (-1, ())
    for p in ps:
        for n in ns:
            for row in endpoint_split(a, p, n):
                ids, child, pc, nc, h, ell, mu, imbalance = row
                edge = max(0, *mu)
                continuation, trace = liberal_path(child)
                candidate = (
                    edge + continuation,
                    ((ids, pc, nc, h, ell, mu, edge),) + trace,
                )
                if candidate[0] > best[0]:
                    best = candidate
    return best


def matrix_from_code(n, code):
    a = [[0] * n for _ in range(n)]
    bit = 0
    for i in range(n):
        for j in range(i + 1, n):
            a[i][j] = a[j][i] = 1 if code >> bit & 1 else -1
            bit += 1
    return tuple(tuple(row) for row in a)


def all_sign(n, sign):
    return tuple(tuple(0 if i == j else sign for j in range(n)) for i in range(n))


def block_parent(core, d, d_sign, cross_row):
    n = len(core)
    a = [[0] * (d + n) for _ in range(d + n)]
    for i in range(d):
        for j in range(i + 1, d):
            a[i][j] = a[j][i] = d_sign
    for i in range(d):
        for j in range(n):
            a[i][d + j] = a[d + j][i] = cross_row[j]
    for i in range(n):
        for j in range(n):
            a[d + i][d + j] = core[i][j]
    return tuple(tuple(row) for row in a)


def root_integral_orientation_capacity(a):
    """Best local capacity using at most one orientation on each shore."""
    _, _, ps, ns = extrema(a)
    best = -1
    witness = None
    for p in ps:
        for n in ns:
            rows = endpoint_split(a, p, n)
            value = sum(max(0, *row[6]) + row[7] for row in rows)
            if value > best:
                best, witness = value, (p, n, rows)
    return best, witness


def check_witnesses():
    k3 = all_sign(3, 1)
    assert extrema(k3)[:2] == (6, 2)
    assert liberal_path(k3)[0] == 6
    negative = extrema(k3)[3][0]
    parent3 = block_parent(k3, 2, 1, negative)
    target3 = (1, 1) + negative
    assert energy(parent3, target3) == 12
    assert max(extrema(parent3)[:2]) == 12
    print("K3: P,N,R,V = 6,2,8,6; K2 parent Q=12 and carries gap 8")

    ledger6 = (
        (0,-1,1,-1,-1,-1),(-1,0,1,-1,1,1),(1,1,0,-1,1,-1),
        (-1,-1,-1,0,1,-1),(-1,1,1,1,0,-1),(-1,1,-1,-1,-1,0),
    )
    assert extrema(ledger6)[:2] == (10, 10)
    assert liberal_path(ledger6)[0] == 12
    print("(10.409): P,N,R,V = 10,10,20,12")

    code36 = matrix_from_code(6, 36)
    assert extrema(code36)[:2] == (14, 22)
    assert liberal_path(code36)[0] == 20
    print("code36: P,N,R,V = 14,22,36,20")
    for row in code36:
        print(" ", row)

    code37 = matrix_from_code(6, 37)
    p37, n37, ps37, ns37 = extrema(code37)
    assert (p37, n37, len(ps37), len(ns37)) == (14, 18, 1, 1)
    integral, witness = root_integral_orientation_capacity(code37)
    mus = tuple(row[6] for row in witness[2])
    assert integral == 24 and sorted(mus) == [(8, 8), (12, 16)]
    assert liberal_path(code37)[0] == 18
    print("code37: P,N,R,V = 14,18,32,18; unique cut mus", mus, "integral local cap", integral)
    for row in code37:
        print(" ", row)

    # Make the negative parent orientation carry p37.  D is negative K5 and
    # every cross row is -p37.  Its target is (1_D,p37), of energy -66.
    inherited = ps37[0]
    parent37 = block_parent(code37, 5, -1, tuple(-z for z in inherited))
    target37 = (1,) * 5 + inherited
    pp, nn, _, _ = extrema(parent37)
    assert (pp, nn, energy(parent37, target37)) == (54, 66, -66)
    assert 18 - (-14) == 32
    print("code37 parent: P,N,Q,target = 54,66,66,-66; carried core gap = 32")

    code8 = matrix_from_code(8, 228544878)
    p8, n8, ps8, ns8 = extrema(code8)
    assert (p8, n8, len(ps8), len(ns8)) == (32, 32, 1, 1)
    value8, trace8 = liberal_path(code8)
    assert value8 == 30
    rows8 = endpoint_split(code8, ps8[0], ns8[0])
    child_data = sorted((r[2], r[3], r[4], r[5], r[6], liberal_path(r[1])[0]) for r in rows8)
    assert child_data == [(4,12,0,16,(20,20),10),(12,4,0,16,(20,20),10)]
    print("code8: P=N=32, unique cut; child data", child_data, "best path", value8, "< R/2=32")
    for row in code8:
        print(" ", row)

    # The same order-eight matrix is the first member of an exact family:
    # C=ss^T-I, sum(s)=sqrt(m), and [[C,J],[J,-C]] for m=4^r.
    s = (1, 1, 1, -1)
    c = tuple(tuple(0 if i == j else s[i] * s[j] for j in range(4)) for i in range(4))
    family4 = tuple(
        tuple(
            (c[i][j] if i < 4 and j < 4 else
             -c[i-4][j-4] if i >= 4 and j >= 4 else 1)
            for j in range(8)
        )
        for i in range(8)
    )
    assert extrema(family4)[:2] == (32, 32)
    assert len(extrema(family4)[2]) == len(extrema(family4)[3]) == 1
    assert liberal_path(family4)[0] == 30
    print("family m=4 base: P=N=32, V=30")


def exhaustive_small():
    expected_strict = {3:(6,8), 4:(10,16), 5:(18,24), 6:(20,36)}
    for n in range(2, 7):
        worst_strict = None
        worst_integral = None
        count = 1 << (n * (n - 1) // 2)
        for code in range(count):
            a = matrix_from_code(n, code)
            p, nn, _, _ = extrema(a)
            r = p + nn
            value, _ = liberal_path(a)
            if p != nn:
                row = (value, r, code)
                if worst_strict is None or value * worst_strict[1] < worst_strict[0] * r:
                    worst_strict = row
            integral, _ = root_integral_orientation_capacity(a)
            row = (integral, r, code)
            if worst_integral is None or integral * worst_integral[1] < worst_integral[0] * r:
                worst_integral = row
        print("exhaustive n", n, "matrices", count, "strict V/R", worst_strict, "integral-root/R", worst_integral)
        if n in expected_strict:
            assert worst_strict[:2] == expected_strict[n]
        if n <= 5:
            assert worst_integral[0] >= worst_integral[1]
        if n == 6:
            assert worst_integral == (24, 32, 37)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--exhaustive-small", action="store_true")
    args = parser.parse_args()
    check_witnesses()
    if args.exhaustive_small:
        exhaustive_small()


if __name__ == "__main__":
    main()
