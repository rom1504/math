#!/usr/bin/env python3
"""Exact checks for Wave 13 joint common-mosaic block replacement."""

from collections import Counter
from functools import lru_cache
from itertools import combinations, product

import numpy as np


@lru_cache(None)
def spins(n):
    return np.array([(1,) + z for z in product((-1, 1), repeat=n - 1)],
                    dtype=np.int16)


def energies(a):
    x = spins(len(a))
    return np.sum((x @ a) * x, axis=1)


def q(a):
    return int(np.max(np.abs(energies(a))))


def gauge_signing(n, mask):
    a = np.zeros((n, n), dtype=np.int16)
    a[0, 1:] = a[1:, 0] = 1
    bit = 0
    for i in range(1, n):
        for j in range(i + 1, n):
            a[i, j] = a[j, i] = 1 if (mask >> bit) & 1 else -1
            bit += 1
    return a


@lru_cache(None)
def gauge_minimizers(n):
    if n == 1:
        return 0, (np.zeros((1, 1), dtype=np.int16),)
    best = n * (n - 1)
    out = []
    for mask in range(1 << ((n - 1) * (n - 2) // 2)):
        a = gauge_signing(n, mask)
        qa = q(a)
        if qa < best:
            best, out = qa, [a]
        elif qa == best:
            out.append(a)
    return best, tuple(out)


@lru_cache(None)
def all_minimizers(n):
    qn, reps = gauge_minimizers(n)
    out = {}
    for a in reps:
        for tail in product((-1, 1), repeat=n - 1):
            d = np.array((1,) + tail, dtype=np.int16)
            b = a * np.outer(d, d)
            out[b.tobytes()] = b
    return qn, tuple(out.values())


def hybrid(a, blocks, replacements):
    b = a.copy()
    for block, g in zip(blocks, replacements):
        b[np.ix_(block, block)] = g
    return b


def best_hybrid(a, blocks):
    choices = [all_minimizers(len(block))[1] for block in blocks]
    best = None
    witnesses = []
    for gs in product(*choices):
        val = q(hybrid(a, blocks, gs))
        if best is None or val < best:
            best, witnesses = val, [gs]
        elif val == best:
            witnesses.append(gs)
    return best, witnesses


A6 = np.array([
    [0, 1, 1, 1, 1, 1],
    [1, 0, -1, -1, 1, 1],
    [1, -1, 0, 1, -1, 1],
    [1, -1, 1, 0, 1, -1],
    [1, 1, -1, 1, 0, -1],
    [1, 1, 1, -1, -1, 0],
], dtype=np.int16)


A7 = np.array([
    [0, 1, 1, 1, 1, 1, 1],
    [1, 0, 1, 1, -1, -1, 1],
    [1, 1, 0, 1, -1, 1, -1],
    [1, 1, 1, 0, 1, -1, -1],
    [1, -1, -1, 1, 0, -1, -1],
    [1, -1, 1, -1, -1, 0, -1],
    [1, 1, -1, -1, -1, -1, 0],
], dtype=np.int16)


A8 = np.array([
    [0, 1, 1, 1, 1, 1, 1, 1],
    [1, 0, 1, -1, 1, 1, -1, -1],
    [1, 1, 0, 1, -1, 1, -1, -1],
    [1, -1, 1, 0, -1, -1, -1, 1],
    [1, 1, -1, -1, 0, -1, 1, -1],
    [1, 1, 1, -1, -1, 0, 1, 1],
    [1, -1, -1, -1, 1, 1, 0, 1],
    [1, -1, -1, 1, -1, 1, 1, 0],
], dtype=np.int16)


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
], dtype=np.int16)


def canonical_bipartitions(n):
    vertices = tuple(range(n))
    for r in range(1, n // 2 + 1):
        for left in combinations(vertices, r):
            if 2 * r == n and 0 not in left:
                continue
            right = tuple(v for v in vertices if v not in left)
            yield (tuple(left), right)


def exhaustive_through_six():
    rows = []
    for n in range(2, 7):
        qn, reps = gauge_minimizers(n)
        worst_best_excess = -1
        positive = 0
        cases = 0
        witness = None
        for a in reps:
            assert q(a) == qn
            for blocks in canonical_bipartitions(n):
                best, _ = best_hybrid(a, blocks)
                excess = best - qn
                cases += 1
                positive += excess > 0
                if excess > worst_best_excess:
                    worst_best_excess = excess
                    witness = (a.copy(), blocks, best)
        rows.append((n, qn, len(reps), cases, positive, worst_best_excess,
                     witness[1]))
    return rows


def selected_order_seven():
    assert q(A7) == 18
    rows = []
    for blocks in canonical_bipartitions(7):
        # Full products are modest for 1+6, 2+5 and 3+4.
        best, witnesses = best_hybrid(A7, blocks)
        rows.append((tuple(map(len, blocks)), blocks[0], best - 18,
                     len(witnesses)))
    return rows


def replacement_profiles(n, blocks):
    """All diagonal replacement energy profiles on gauge-fixed global spins."""
    x = spins(n)
    local = []
    for block in blocks:
        xb = x[:, block]
        profiles = []
        for g in all_minimizers(len(block))[1]:
            profiles.append(np.sum((xb @ g) * xb, axis=1))
        local.append(np.array(profiles, dtype=np.int16))
    sums = (local[0][:, None, :] + local[1][None, :, :]).reshape(
        -1, len(x))
    return np.unique(sums, axis=0)


def exhaustive_order_seven():
    """Fixed canonical shores suffice because all labelled gauges are scanned."""
    qn, reps = gauge_minimizers(7)
    assert qn == 18 and len(reps) == 3240
    out = []
    x = spins(7)
    for r in (1, 2, 3):
        blocks = (tuple(range(r)), tuple(range(r, 7)))
        profiles = replacement_profiles(7, blocks)
        positive = 0
        max_excess = -1
        excess_counts = {}
        for a in reps:
            cross = a.copy()
            for block in blocks:
                cross[np.ix_(block, block)] = 0
            ce = np.sum((x @ cross) * x, axis=1)
            vals = np.max(np.abs(profiles + ce[None, :]), axis=1)
            excess = int(np.min(vals)) - qn
            excess_counts[excess] = excess_counts.get(excess, 0) + 1
            positive += excess > 0
            max_excess = max(max_excess, excess)
        out.append((tuple(map(len, blocks)), len(profiles), positive,
                    max_excess, tuple(sorted(excess_counts.items()))))
    return out


def selected_order_eight():
    assert q(A8) == 20
    x = spins(8)
    out = []
    for r in (2, 3, 4):
        blocks = (tuple(range(r)), tuple(range(r, 8)))
        profiles = replacement_profiles(8, blocks)
        cross = A8.copy()
        for block in blocks:
            cross[np.ix_(block, block)] = 0
        ce = np.sum((x @ cross) * x, axis=1)
        vals = np.max(np.abs(profiles + ce[None, :]), axis=1)
        best = int(np.min(vals))
        out.append((tuple(map(len, blocks)), len(profiles), best - 20,
                    int(np.sum(vals == best)), tuple(sorted(set(map(int, vals))))))
    return out


def selected_order_nine():
    assert q(A9) == 24
    x = spins(9)
    out = []
    for r in (3, 4):
        blocks = (tuple(range(r)), tuple(range(r, 9)))
        profiles = replacement_profiles(9, blocks)
        cross = A9.copy()
        for block in blocks:
            cross[np.ix_(block, block)] = 0
        ce = np.sum((x @ cross) * x, axis=1)
        vals = np.max(np.abs(profiles + ce[None, :]), axis=1)
        best = int(np.min(vals))
        detail = None
        if r == 3:
            dqs = tuple(q(A9[np.ix_(block, block)]) for block in blocks)
            detail = (dqs, tuple(all_minimizers(len(block))[0] for block in blocks),
                      q(cross), int(np.sum(ce == q(cross))),
                      int(np.sum(ce == -q(cross))))
        out.append((tuple(map(len, blocks)), len(profiles), best - 24,
                    int(np.sum(vals == best)), tuple(sorted(set(map(int, vals)))),
                    detail))
    return out


def order_nine_obstruction_witness():
    blocks = (tuple(range(3)), tuple(range(3, 9)))
    choices = [all_minimizers(len(block))[1] for block in blocks]
    best, witness = best_hybrid(A9, blocks)
    assert best == 28
    gs = witness[0]
    b = hybrid(A9, blocks, gs)
    es = energies(b)
    return gs, b, tuple((int(v), int(np.sum(es == v)))
                        for v in sorted(set(map(int, es))))


def order_nine_six_state_certificate():
    """A compact simultaneous-state certificate that Q<=24 is impossible."""
    xa = np.array((1, -1, 1), dtype=np.int16)
    xb = np.array((1, 1, 1), dtype=np.int16)
    ya = np.array((1, 1, 1, 1, 1, -1), dtype=np.int16)
    yb = np.array((1, 1, 1, 1, -1, -1), dtype=np.int16)

    def e(a, z):
        return int(z @ a @ z)

    p3 = Counter((e(g, xa), e(g, xb)) for g in all_minimizers(3)[1])
    p6 = Counter((e(g, ya), e(g, yb)) for g in all_minimizers(6)[1])
    expected3 = Counter({(-6, 2): 1, (-2, -2): 2, (-2, 6): 1,
                         (2, -6): 1, (2, 2): 2, (6, -2): 1})
    expected6 = Counter({(-10, -6): 60, (-10, 10): 12,
                         (-6, -10): 60, (-6, 6): 60,
                         (6, -6): 60, (6, 10): 60,
                         (10, -10): 12, (10, 6): 60})
    assert p3 == expected3 and p6 == expected6

    # The six global states are the two signs of three projective pairs.
    # Q<=24 would force a+u=0, |a+v|<=4 and |b+v|<=4.
    feasible = []
    for (a, b), n3 in p3.items():
        for (u, v), n6 in p6.items():
            if a + u == 0 and abs(a + v) <= 4 and abs(b + v) <= 4:
                feasible.append(((a, b), (u, v), n3 * n6))
    assert feasible == []

    # Independently verify the six states cover all 8*384 choices.
    states = []
    for x0, y0 in ((xa, ya), (xa, yb), (xb, yb)):
        states.extend((np.concatenate((x0, y0)),
                       np.concatenate((x0, -y0))))
    choices3 = all_minimizers(3)[1]
    choices6 = all_minimizers(6)[1]
    assert (len(choices3), len(choices6)) == (8, 384)
    for g3, g6 in product(choices3, choices6):
        h = hybrid(A9, (tuple(range(3)), tuple(range(3, 9))), (g3, g6))
        assert any(abs(int(z @ h @ z)) > 24 for z in states)

    # No single state is an obstruction: if the minimizers may depend on z,
    # every state can be driven to absolute energy at most eight.
    x = spins(9)
    blocks = (tuple(range(3)), tuple(range(3, 9)))
    profiles = replacement_profiles(9, blocks)
    cross = A9.copy()
    for block in blocks:
        cross[np.ix_(block, block)] = 0
    ce = np.sum((x @ cross) * x, axis=1)
    assert int(np.max(np.min(np.abs(profiles + ce[None, :]), axis=0))) == 8
    return p3, p6


def j_value(c, a):
    y = spins(len(a))
    exposure = 2 * np.sum(np.abs(y @ c.T), axis=1)
    return int(np.max(exposure + np.abs(energies(a))))


def check_order_nine_local_vs_global():
    left, right = tuple(range(3)), tuple(range(3, 9))
    c = A9[np.ix_(left, right)]
    d3 = A9[np.ix_(left, left)]
    d6 = A9[np.ix_(right, right)]
    tail_vals = [j_value(c, g) for g in all_minimizers(6)[1]]
    head_vals = [j_value(c.T, g) for g in all_minimizers(3)[1]]
    assert (q(d6), j_value(c, d6), min(tail_vals)) == (14, 30, 30)
    assert (q(d3), j_value(c.T, d3), min(head_vals)) == (6, 30, 26)
    # Common-mosaic responses: rise = internal excess + hybrid excess.
    assert 24 - q(d3) - q(d6) == 4
    assert 28 - all_minimizers(3)[0] - all_minimizers(6)[0] == 12


def check_deterministic_bound():
    # Exhaust every n<=6 minimizing case and every candidate replacement in
    # the returned best witness, checking both block subadditivity and 2n(s-1).
    for n in range(2, 7):
        qn, reps = gauge_minimizers(n)
        for a in reps:
            for blocks in canonical_bipartitions(n):
                best, witnesses = best_hybrid(a, blocks)
                gs = witnesses[0]
                hs = []
                for block, g in zip(blocks, gs):
                    d = a[np.ix_(block, block)]
                    hs.append(g - d)
                h = np.zeros_like(a)
                for block, hb in zip(blocks, hs):
                    h[np.ix_(block, block)] = hb
                assert q(h) <= sum(q(hb) for hb in hs)
                assert best - qn <= q(h)
                s = max(map(len, blocks))
                assert q(h) <= 2 * n * (s - 1)


if __name__ == "__main__":
    print("all minimizing gauge representatives and bipartitions, n<=6")
    for row in exhaustive_through_six():
        print(row)
    rows7 = selected_order_seven()
    assert all(row[2] == 0 for row in rows7)
    print("selected A7: all 63 bipartitions have zero best excess")
    print("all 3240 order-7 minimizing gauges, one canonical split per size")
    for row in exhaustive_order_seven():
        print(row)
    print("selected A8 canonical 2+6, 3+5, 4+4 splits")
    for row in selected_order_eight():
        print(row)
    print("selected A9 canonical 3+6 and 4+5 splits")
    for row in selected_order_nine():
        print(row)
    gs9, b9, hist9 = order_nine_obstruction_witness()
    print("A9 3+6 achieving minimizers:")
    print(gs9[0])
    print(gs9[1])
    print("hybrid energy histogram:", hist9)
    p3, p6 = order_nine_six_state_certificate()
    print("A9 six-state profile tables:", sorted(p3.items()), sorted(p6.items()))
    print("A9 six-state simultaneous obstruction: PASS")
    check_order_nine_local_vs_global()
    print("A9 local-versus-global separation: PASS")
    check_deterministic_bound()
    print("deterministic bound checks: PASS")
