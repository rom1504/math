#!/usr/bin/env python3
"""Independent exact audit of the Wave 13 A8 residual obstruction."""

from fractions import Fraction
from itertools import product

import numpy as np


A8 = np.array([
    [0, 1, 1, 1, 1, 1, 1, 1],
    [1, 0, 1, -1, 1, 1, -1, -1],
    [1, 1, 0, 1, -1, 1, -1, -1],
    [1, -1, 1, 0, -1, -1, -1, 1],
    [1, 1, -1, -1, 0, -1, 1, -1],
    [1, 1, 1, -1, -1, 0, 1, 1],
    [1, -1, -1, -1, 1, 1, 0, 1],
    [1, -1, -1, 1, -1, 1, 1, 0],
], dtype=np.int64)


def projective_spins(n):
    for tail in product((-1, 1), repeat=n - 1):
        yield np.array((1,) + tail, dtype=np.int64)


def energy(a, x):
    return int(x @ a @ x)


def pnq(a):
    values = [(energy(a, x), x) for x in projective_spins(len(a))]
    p = max(v for v, _ in values)
    n = -min(v for v, _ in values)
    return p, n, max(p, n), values


def residual_data(a, p, n, xidx):
    xidx = tuple(xidx)
    didx = tuple(i for i in range(len(a)) if i not in xidx)
    x = p[list(xidx)]
    d = p[list(didx)]
    ax = a[np.ix_(xidx, xidx)]
    ad = a[np.ix_(didx, didx)]
    b = a[np.ix_(didx, xidx)]
    px, nx, qx, _ = pnq(ax)
    hx = energy(ax, x)
    lx = int(np.sum(np.abs(b @ x)))
    parent_q = pnq(a)[2]
    decrement = parent_q - qx
    capacities = []
    residuals = []
    for sigma in (1, -1):
        raw = 2 * lx - (qx - sigma * hx)
        cap = max(0, raw)
        residual = max(0, cap - decrement)
        direct = max(0, 2 * lx + sigma * hx - parent_q)
        assert residual == direct
        capacities.append(cap)
        residuals.append(residual)
        if residual:
            by = b @ x
            u = np.where(sigma * by >= 0, 1, -1)
            assert sigma * energy(ad, u) <= -residual
    return {
        "X": xidx,
        "D": didx,
        "P_X": px,
        "N_X": nx,
        "Q_X": qx,
        "h_X": hx,
        "L_X": lx,
        "decrement": decrement,
        "capacities": tuple(capacities),
        "residuals": tuple(residuals),
        "D_ground_energy": energy(ad, d),
    }


def exact_event_probability(row_fields, deleted):
    """Probability of exactly `deleted` when H is that shore in (10.330)."""
    out = Fraction(1)
    for i in deleted:
        out *= Fraction(int(row_fields[i]), len(row_fields) - 1)
    return out


def main():
    pval, nval, qval, vals = pnq(A8)
    positives = [x for v, x in vals if v == pval]
    negatives = [x for v, x in vals if v == -nval]
    assert (pval, nval, qval, len(positives), len(negatives)) == (20, 20, 20, 4, 4)

    summaries = []
    admissible = []
    for p in positives:
        switched_rows = p * (A8 @ p)
        assert np.all(switched_rows >= 0)
        assert int(np.sum(switched_rows)) == 20
        for n in negatives:
            labels = p * n
            shores = (tuple(np.flatnonzero(labels == 1)), tuple(np.flatnonzero(labels == -1)))
            assert tuple(sorted(map(len, shores))) == (4, 4)
            pair = []
            for xidx in shores:
                datum = residual_data(A8, p, n, xidx)
                assert datum["P_X"] == datum["N_X"] == datum["Q_X"] == 8
                assert datum["h_X"] == 0
                assert datum["L_X"] == 10
                assert datum["decrement"] == 12
                assert datum["capacities"] == (12, 12)
                assert datum["residuals"] == (0, 0)
                pair.append(datum)

                deleted = datum["D"]
                prob = exact_event_probability(switched_rows, deleted)
                if prob:
                    admissible.append((tuple(map(int, switched_rows)), deleted, datum["X"], prob))
            summaries.append((tuple(map(int, switched_rows)), shores, pair))

    assert len(summaries) == 16
    row_profiles = sorted(set(row for row, _, _ in summaries))
    print("P,N,Q and projective endpoint counts:", pval, nval, qval, len(positives), len(negatives))
    print("positive-endpoint switched row profiles:")
    for row in row_profiles:
        print(" ", row)
    print("endpoint pairs audited:", len(summaries))
    print("all shores: size 4, P=N=Q=8, h=0, L=10, capacity=(12,12), residual=(0,0)")
    print("field-proportional exact-shore outcomes with positive probability:", len(admissible))
    for item in admissible[:12]:
        print(" ", item)
    print("centered deterministic/outcome gap: (20-5*sqrt(2))-12 = 8-5*sqrt(2) > 0")


if __name__ == "__main__":
    main()
