#!/usr/bin/env python3
"""Exact/numerical checks for the Wave 42 mixed-pressure audit."""

from __future__ import annotations

import itertools
import math

import numpy as np

from envelope_block_cover_r27 import A6, A8, A9


def spins(n):
    return [np.array((1,) + t, dtype=np.int64) for t in itertools.product((-1, 1), repeat=n - 1)]


def qnorm(a, ys=None):
    if ys is None:
        ys = spins(len(a))
    return max(abs(int(y @ a @ y)) for y in ys)


def selector_qs(a, m):
    ys = spins(m)
    return {
        s: qnorm(a[np.ix_(s, s)], ys)
        for s in itertools.combinations(range(len(a)), m)
    }


def check_q_self_bounding():
    for a in (A6, A8, A9):
        n = len(a)
        for m in range(2, n + 1):
            qs = selector_qs(a, m)
            qprev = selector_qs(a, m - 1) if m > 1 else {}
            edge_square_sum = 0
            edge_count = 0
            for s, q in qs.items():
                # Pick a ground and its sector, then verify local stability and
                # the deletion inequalities exactly.
                b = a[np.ix_(s, s)]
                ground = next(y for y in spins(m) if abs(int(y @ b @ y)) == q)
                raw = int(ground @ b @ ground)
                eps = 1 if raw >= 0 else -1
                local = [eps * int(ground[j] * (b @ ground)[j]) for j in range(m)]
                assert min(local) >= 0 and sum(local) == q
                deletions = []
                for pos, i in enumerate(s):
                    core = tuple(j for j in s if j != i)
                    delta = q - qprev[core]
                    deletions.append(delta)
                    assert 0 <= delta <= 2 * local[pos]
                assert sum(deletions) <= 2 * q
                assert sum(d * d for d in deletions) <= 4 * (m - 1) * q

                if m < n:
                    for i in s:
                        for j in range(n):
                            if j in s:
                                continue
                            t = tuple(sorted((set(s) - {i}) | {j}))
                            edge_square_sum += (qs[t] - q) ** 2
                            edge_count += 1

            if m < n:
                vals = list(qs.values())
                mean = sum(vals) / len(vals)
                var = sum((v - mean) ** 2 for v in vals) / len(vals)
                edge_mean = edge_square_sum / edge_count
                poincare_rhs = m * (n - m) * edge_mean / (2 * n)
                assert var <= poincare_rhs + 1e-12
                assert edge_mean <= 8 * mean + 1e-12


def check_deficit_non_self_bound():
    a = A8
    z = spins(8)[0]
    sigma = -1
    s = (0, 1, 3, 4, 5, 6)
    qs = selector_qs(a, 6)
    q5 = selector_qs(a, 5)
    b = a[np.ix_(s, s)]
    c = int(z[list(s)] @ b @ z[list(s)])
    d = qs[s] - sigma * c
    increments = []
    for i in s:
        core = tuple(j for j in s if j != i)
        bc = a[np.ix_(core, core)]
        cc = int(z[list(core)] @ bc @ z[list(core)])
        dc = q5[core] - sigma * cc
        increments.append(d - dc)
    assert (qs[s], c, d) == (18, -14, 4)
    assert increments == [-4, 4, -4, 4, 4, 4]
    assert sum(max(x, 0) for x in increments) == 4 * d


def check_pressure_identities():
    a = A8
    n, m = 8, 4
    z = spins(n)[0]
    sigma = -1
    qs = selector_qs(a, m)
    ds = []
    for s, q in qs.items():
        b = a[np.ix_(s, s)]
        c = int(z[list(s)] @ b @ z[list(s)])
        d = q - sigma * c
        assert d >= 0
        ds.append(d)

    theta, h = 0.2, 10.0
    w = np.exp(theta * (h - np.array(ds, dtype=float)))
    k1 = math.log(float(np.mean(w)))
    k2 = math.log(float(np.mean(w * w)))
    gap = k2 - 2 * k1
    density = w / float(np.mean(w))
    renyi = math.log(float(np.mean(density * density)))
    assert abs(gap - renyi) < 1e-13
    assert gap <= theta * h - k1 + 1e-13

    favorable = sum(d <= h for d in ds) / len(ds)
    direct = (math.exp(k1) - 1.0) / (math.exp(theta * h) - 1.0)
    assert k1 > 0 and favorable + 1e-13 >= direct

    # Exact edge/core port formula for every adjacent selector pair.
    q3 = selector_qs(a, 3)
    for s in qs:
        for aa in s:
            core = tuple(i for i in s if i != aa)
            for bb in range(n):
                if bb in s:
                    continue
                t = tuple(sorted(core + (bb,)))
                cs = int(z[list(s)] @ a[np.ix_(s, s)] @ z[list(s)])
                ct = int(z[list(t)] @ a[np.ix_(t, t)] @ z[list(t)])
                cc = int(z[list(core)] @ a[np.ix_(core, core)] @ z[list(core)])
                ds0 = qs[s] - sigma * cs
                dt0 = qs[t] - sigma * ct
                ra = (qs[s] - q3[core]) - sigma * (cs - cc)
                rb = (qs[t] - q3[core]) - sigma * (ct - cc)
                assert dt0 - ds0 == rb - ra


if __name__ == "__main__":
    check_q_self_bounding()
    check_deficit_non_self_bound()
    check_pressure_identities()
    print("PASS principal-norm deletion self-bound and Johnson Poincare audit")
    print("PASS mixed-pressure Renyi/cap/direct-tail identities")
    print("PASS exact A8 deficit non-self-bound and selector port formula")
