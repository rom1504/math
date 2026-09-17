#!/usr/bin/env python3
"""Exact one-deletion census for the Wave-17 selected-child route."""

from collections import Counter
from fractions import Fraction
from itertools import combinations, product


A9 = (
    (0, 1, -1, 1, 1, 1, 1, -1, -1),
    (1, 0, 1, -1, -1, -1, 1, -1, -1),
    (-1, 1, 0, 1, 1, 1, 1, 1, -1),
    (1, -1, 1, 0, 1, 1, 1, -1, 1),
    (1, -1, 1, 1, 0, 1, -1, 1, -1),
    (1, -1, 1, 1, 1, 0, -1, -1, -1),
    (1, 1, 1, 1, -1, -1, 0, 1, 1),
    (-1, -1, 1, -1, 1, -1, 1, 0, -1),
    (-1, -1, -1, 1, -1, -1, 1, -1, 0),
)

QMIN = {1: 0, 2: 2, 3: 6, 4: 8, 5: 8, 6: 10, 7: 18, 8: 20, 9: 24}


def data(vertices):
    vertices = tuple(vertices)
    n = len(vertices)
    if n <= 1:
        return 0, 1, Counter({0: 1}), ((1,),)
    vals = []
    spins = []
    for tail in product((-1, 1), repeat=n - 1):
        x = (1,) + tail
        e = sum(
            2 * A9[vertices[i]][vertices[j]] * x[i] * x[j]
            for i in range(n) for j in range(i + 1, n)
        )
        vals.append(e)
        spins.append(x)
    q = max(map(abs, vals))
    grounds = tuple(x for x, e in zip(spins, vals) if abs(e) == q)
    deficits = Counter(q - abs(e) for e in vals)
    return q, len(grounds), deficits, grounds


def deletion_record(vertices):
    q, ng, deficits, grounds = data(vertices)
    children = []
    for v in vertices:
        child = tuple(x for x in vertices if x != v)
        qc, ngc, dc, gc = data(child)
        children.append((v, qc, qc - QMIN[len(child)], ngc, dc))
    return q, q - QMIN[len(vertices)], ng, deficits, children


def main():
    allv = tuple(range(9))
    chain = [allv, tuple(v for v in allv if v != 7)]
    # Take the lexicographically first order-seven child of the flat B.
    chain.append(chain[-1][1:])
    for vertices in chain:
        q, eps, ng, deficits, children = deletion_record(vertices)
        print("V=", vertices, "r=", len(vertices), "Q=", q,
              "eps=", eps, "grounds(projective)=", ng,
              "deficits=", sorted(deficits.items()))
        print(" children (v,Q,eps,#grounds):",
              [(v, qc, ec, ngc) for v, qc, ec, ngc, _ in children])

    # Every order-r restriction: joint histogram of parent excess, grounds,
    # minimum child excess and child ground counts.
    for r in range(3, 10):
        hist = Counter()
        for vertices in combinations(allv, r):
            q, eps, ng, _, children = deletion_record(vertices)
            min_eps = min(c[2] for c in children)
            max_dec = max(q - c[1] for c in children)
            hist[(eps, ng, min_eps, max_dec)] += 1
        print("order", r, "state histogram", sorted(hist.items()))

    # Exact normalized increments on the named flat chain.
    B = chain[1]
    qB, epsB, *_ = deletion_record(B)
    for v, qc, ec, ngc, _ in deletion_record(B)[4]:
        inc2 = Fraction(ec * ec, 7**3) - Fraction(epsB * epsB, 8**3)
        print("flat 8->7 v", v, "e^2 increment", inc2)


if __name__ == "__main__":
    main()
