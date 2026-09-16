#!/usr/bin/env python3
"""Exact audits for two-directed-layer endpoint partition trees.

At every endpoint-generated split U=S union T, select the better of the two
Q-benchmarked orientations toward each shore and recurse into both shores.
The resulting value is deliberately liberal: endpoint ties are optimized
independently at every node.  Energies use doubled normalization.
"""

from functools import lru_cache
from itertools import product


def energy(a, x):
    return sum(a[i][j] * x[i] * x[j]
               for i in range(len(a)) for j in range(len(a)))


def induced(a, ids):
    return tuple(tuple(a[i][j] for j in ids) for i in ids)


@lru_cache(None)
def extrema(a):
    if not a:
        return 0, 0, ((),), ((),)
    xs = (bits + (1,) for bits in product((-1, 1), repeat=len(a) - 1))
    rows = tuple((energy(a, x), x) for x in xs)
    hi = max(v for v, _ in rows)
    lo = min(v for v, _ in rows)
    return (hi, -lo, tuple(x for v, x in rows if v == hi),
            tuple(x for v, x in rows if v == lo))


def pair_rows(a, p, n):
    sides = [tuple(i for i in range(len(a)) if p[i] * n[i] == bit)
             for bit in (-1, 1)]
    if not all(sides):
        return ()
    out = []
    for xids in sides:
        yids = tuple(i for i in range(len(a)) if i not in xids)
        ax = induced(a, xids)
        px = tuple(p[i] for i in xids)
        h = energy(ax, px)
        pp, nn, _, _ = extrema(ax)
        q = max(pp, nn)
        fields = [sum(a[j][i] * p[i] for i in xids) for j in yids]
        ell = sum(abs(v) for v in fields)
        mus = (max(0, 2 * ell - (q - h)),
               max(0, 2 * ell - (q + h)))
        out.append((xids, pp, nn, h, ell, mus, ax))
    return tuple(out)


@lru_cache(None)
def branch_value(a):
    """Maximum full branching capacity and a witness tree."""
    if len(a) <= 1:
        return 0, None
    pp, nn, ps, ns = extrema(a)
    best = (-1, None)
    for p in ps:
        for n in ns:
            rows = pair_rows(a, p, n)
            if len(rows) != 2:
                continue
            child = [branch_value(row[-1]) for row in rows]
            here = sum(max(row[5]) for row in rows)
            value = here + sum(v for v, _ in child)
            if value > best[0]:
                best = (value, (p, n, rows, child))
    return best


def summarize(name, a):
    pp, nn, _, _ = extrema(a)
    value, witness = branch_value(a)
    print(name, "n", len(a), "P", pp, "N", nn, "R", pp + nn,
          "Q", max(pp, nn), "T", value,
          "T/Q", value / max(pp, nn) if max(pp, nn) else 0,
          "T/R", value / (pp + nn) if pp + nn else 0)
    if witness:
        for row, child in zip(witness[2], witness[3]):
            print(" child", row[:-1], "Tchild", child[0])


def clique(n, sign):
    return tuple(tuple(0 if i == j else sign for j in range(n))
                 for i in range(n))


A6 = (
    (0, -1, 1, -1, -1, -1),
    (-1, 0, 1, -1, 1, 1),
    (1, 1, 0, -1, 1, -1),
    (-1, -1, -1, 0, 1, -1),
    (-1, 1, 1, 1, 0, -1),
    (-1, 1, -1, -1, -1, 0),
)

A12 = (
    (0,-1,-1,-1,-1,-1,-1,1,-1,1,-1,-1),
    (-1,0,-1,-1,-1,-1,1,-1,1,-1,-1,-1),
    (-1,-1,0,1,-1,1,1,-1,1,1,1,-1),
    (-1,-1,1,0,-1,1,1,1,1,1,1,1),
    (-1,-1,-1,-1,0,1,1,1,1,-1,1,1),
    (-1,-1,1,1,1,0,1,-1,-1,-1,-1,-1),
    (-1,1,1,1,1,1,0,-1,1,-1,-1,-1),
    (1,-1,-1,1,1,-1,-1,0,-1,1,1,1),
    (-1,1,1,1,1,-1,1,-1,0,1,1,1),
    (1,-1,1,1,-1,-1,-1,1,1,0,1,-1),
    (-1,-1,1,1,1,-1,-1,1,1,1,0,-1),
    (-1,-1,-1,1,1,-1,-1,1,1,-1,-1,0),
)


if __name__ == "__main__":
    for n in range(2, 13):
        summarize(f"K+{n}", clique(n, 1))
        summarize(f"K-{n}", clique(n, -1))
    summarize("ledger-A6", A6)
    summarize("reset-A12", A12)
