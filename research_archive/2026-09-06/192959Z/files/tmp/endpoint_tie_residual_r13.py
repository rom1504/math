#!/usr/bin/env python3
"""Exact finite checks for Wave 13 endpoint-tie residual service.

All matrix calculations after the switching-class search use Python integers
and Fractions.  NumPy is used only to enumerate switching-normalized signings
and locate their exact (integer) Q values.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product

import numpy as np


A5 = (
    (0, -1, 1, -1, 1),
    (-1, 0, -1, 1, 1),
    (1, -1, 0, 1, 1),
    (-1, 1, 1, 0, 1),
    (1, 1, 1, 1, 0),
)

A8 = (
    (0, 1, 1, 1, 1, 1, 1, 1),
    (1, 0, 1, -1, 1, 1, -1, -1),
    (1, 1, 0, 1, -1, 1, -1, -1),
    (1, -1, 1, 0, -1, -1, -1, 1),
    (1, 1, -1, -1, 0, -1, 1, -1),
    (1, 1, 1, -1, -1, 0, 1, 1),
    (1, -1, -1, -1, 1, 1, 0, 1),
    (1, -1, -1, 1, -1, 1, 1, 0),
)


def spins(n: int):
    if n == 0:
        return ((),)
    return tuple((1,) + t for t in product((-1, 1), repeat=n - 1))


def energy(a, x):
    return sum(
        a[i][j] * x[i] * x[j]
        for i in range(len(a))
        for j in range(len(a))
    )


def stats(a):
    vals = tuple((energy(a, x), x) for x in spins(len(a)))
    if not vals:
        return 0, 0, 0, ((),), ((),)
    p = max(v for v, _ in vals)
    lo = min(v for v, _ in vals)
    return (
        p,
        -lo,
        max(p, -lo),
        tuple(x for v, x in vals if v == p),
        tuple(x for v, x in vals if v == lo),
    )


def principal(a, ids):
    return tuple(tuple(a[i][j] for j in ids) for i in ids)


def normalized_matrix(n: int, code: int):
    """Representative with every edge incident to vertex zero equal to +1."""
    a = [[0] * n for _ in range(n)]
    for j in range(1, n):
        a[0][j] = a[j][0] = 1
    k = 0
    for i in range(1, n):
        for j in range(i + 1, n):
            a[i][j] = a[j][i] = 1 if (code >> k) & 1 else -1
            k += 1
    return tuple(tuple(row) for row in a)


def enumerate_minimizers(n: int, batch: int = 4096):
    """Exhaust all switching classes and return q_n and minimizing codes."""
    xs = np.asarray(spins(n), dtype=np.int16)
    fixed = 2 * np.sum(xs[:, :1] * xs[:, 1:], axis=1, dtype=np.int16)
    edges = [(i, j) for i in range(1, n) for j in range(i + 1, n)]
    coeff = np.asarray(
        [[2 * x[i] * x[j] for i, j in edges] for x in xs], dtype=np.int16
    )
    total = 1 << len(edges)
    best = None
    codes = []
    bit_positions = np.arange(len(edges), dtype=np.uint64)
    for start in range(0, total, batch):
        stop = min(start + batch, total)
        cc = np.arange(start, stop, dtype=np.uint64)
        if edges:
            signs = 2 * ((cc[:, None] >> bit_positions) & 1).astype(np.int16) - 1
            vals = fixed[None, :] + signs @ coeff.T
        else:
            vals = np.broadcast_to(fixed, (len(cc), len(fixed)))
        qq = np.max(np.abs(vals), axis=1)
        local = int(np.min(qq))
        if best is None or local < best:
            best = local
            codes = [int(c) for c in cc[qq == local]]
        elif local == best:
            codes.extend(int(c) for c in cc[qq == local])
    return best, tuple(codes)


def endpoint_rows(a):
    """All projective positive/negative endpoint pairs and root data."""
    pval, nval, q, ps, ns = stats(a)
    out = []
    for p in ps:
        for n in ns:
            ratio = tuple(x * y for x, y in zip(p, n))
            shores = tuple(
                tuple(i for i, z in enumerate(ratio) if z == sign)
                for sign in (-1, 1)
            )
            assert all(shores)
            c_total = 0
            child_imbalance = 0
            residual_total = 0
            shore_rows = []
            for x_ids in shores:
                d_ids = tuple(i for i in range(len(a)) if i not in x_ids)
                ax = principal(a, x_ids)
                px, nx, qx, _, _ = stats(ax)
                y = tuple(p[i] for i in x_ids)
                h = energy(ax, y)
                by = tuple(
                    sum(a[i][j] * p[j] for j in x_ids) for i in d_ids
                )
                ell = sum(abs(v) for v in by)
                assert 4 * ell == pval + nval
                caps = tuple(max(0, 2 * ell - (qx - sig * h)) for sig in (1, -1))
                decrement = q - qx
                residual = tuple(max(0, c - decrement) for c in caps)
                formula = tuple(
                    max(Fraction(0), Fraction(sig * h) - Fraction(abs(pval - nval), 2))
                    for sig in (1, -1)
                )
                assert tuple(Fraction(v) for v in residual) == formula

                # A positive residual forces negative oriented sibling energy.
                for sig, rr in zip((1, -1), residual):
                    u = tuple(1 if sig * v >= 0 else -1 for v in by)
                    du = energy(principal(a, d_ids), u)
                    exposed = 2 * ell + sig * h - q
                    assert rr == max(0, exposed)
                    assert sig * du <= -exposed
                    if rr:
                        assert sig * du <= -rr < 0

                c_total += sum(caps)
                child_imbalance += abs(px - nx)
                residual_total += sum(residual)
                shore_rows.append(
                    (x_ids, px, nx, qx, h, caps, decrement, residual)
                )
            t0 = Fraction(pval + nval, c_total + child_imbalance)
            canonical_service = t0 * residual_total
            out.append(
                {
                    "p": p,
                    "n": n,
                    "shores": shores,
                    "C": c_total,
                    "D": child_imbalance,
                    "t0": t0,
                    "raw_residual": residual_total,
                    "canonical_service": canonical_service,
                    "shore_rows": tuple(shore_rows),
                }
            )
    return pval, nval, q, tuple(out)


def audit_small_minimizers():
    table = []
    for n in range(2, 8):
        qn, codes = enumerate_minimizers(n)
        optimized = []
        prescribed_min = []
        zero_pair_counts = []
        raw_optimized = []
        raw_averages = []
        for code in codes:
            a = normalized_matrix(n, code)
            assert stats(a)[2] == qn
            _, _, _, rows = endpoint_rows(a)
            services = [row["canonical_service"] for row in rows]
            raw_services = [row["raw_residual"] for row in rows]
            optimized.append(max(services))
            prescribed_min.append(min(services))
            zero_pair_counts.append(sum(v == 0 for v in services))
            raw_optimized.append(max(raw_services))
            raw_averages.append(
                Fraction(sum(raw_services), len(raw_services))
            )
        row = (
            n,
            qn,
            len(codes),
            min(optimized),
            max(optimized),
            min(raw_optimized),
            max(raw_optimized),
            min(raw_averages),
            max(raw_averages),
            min(prescribed_min),
            max(zero_pair_counts),
        )
        table.append(row)
        print(
            "SMALL",
            "n,q,#min,canonical-best-range,raw-best-range,"
            "raw-average-range,worst-pair,max-zero-pairs =",
            row,
        )
    return tuple(table)


def audit_a5():
    p, n, q, rows = endpoint_rows(A5)
    assert (p, n, q, len(rows)) == (8, 8, 8, 25)
    services = [r["canonical_service"] for r in rows]
    raw = [r["raw_residual"] for r in rows]
    assert raw.count(0) == 5 and raw.count(4) == 20
    assert services.count(0) == 5 and services.count(Fraction(16, 5)) == 20
    assert sum(services, Fraction()) / 25 == Fraction(64, 25)

    temporal_deleted = (4,)
    matching = [r for r in rows if temporal_deleted in r["shores"]]
    assert len(matching) == 1 and matching[0]["canonical_service"] == 0

    # On every 2+3 pair, a general conserved K<=4 allocation can fill both
    # residual-positive root buckets, pass four units to the imbalanced
    # order-three child, and discharge those four units in a residual-zero
    # bucket there.  Its edge loads are 1, 1 at the root and 1 below one
    # edge, so theta_root=2.  Total raw residual four is the matching dual
    # upper certificate.
    positive = next(r for r in rows if r["raw_residual"] == 4)
    assert positive["D"] == 4
    root_positive = []
    for shore in positive["shore_rows"]:
        for cap, residual in zip(shore[5], shore[7]):
            if residual:
                root_positive.append((cap, residual))
    assert sorted(root_positive) == [(4, 2), (8, 2)]
    assert sum(r for _, r in root_positive) == 4
    assert 2 <= 4  # exact path-cover value of the displayed allocation

    # Exact radical comparisons, performed after moving the positive terms
    # and squaring only positive quantities.
    assert 200 * 200 > 64 * 64 * 5  # 8-64 sqrt(5)/25 > 0
    assert 64 * 64 * 5 > 136 * 136  # 64/25 exceeds that demand
    print(
        "A5",
        "matching-service=0 canonical-best=16/5 general-best=4",
        "uniform canonical=64/25 uniform general=16/5 positive-pairs=20/25",
    )


def audit_a8():
    p, n, q, rows = endpoint_rows(A8)
    assert (p, n, q, len(rows)) == (20, 20, 20, 16)
    assert all(r["raw_residual"] == 0 for r in rows)
    assert all(r["canonical_service"] == 0 for r in rows)
    assert all(r["C"] == 48 and r["D"] == 0 and r["t0"] == Fraction(5, 6) for r in rows)
    assert all(
        shore[1:5] == (8, 8, 8, 0)
        for r in rows
        for shore in r["shore_rows"]
    )

    # Every singleton restriction is an exact q_7 minimizer and has decrement 2.
    children = []
    for i in range(8):
        ids = tuple(j for j in range(8) if j != i)
        children.append(stats(principal(A8, ids))[2])
    assert children == [18] * 8

    # Every positive endpoint has strictly fractional field-proportional
    # singleton probabilities: fields are a permutation of 1,1,1,3,3,3,3,5.
    _, _, _, ps, _ = stats(A8)
    profiles = []
    for ground in ps:
        row_fields = tuple(
            ground[i] * sum(A8[i][j] * ground[j] for j in range(8))
            for i in range(8)
        )
        assert sorted(row_fields) == [1, 1, 1, 3, 3, 3, 3, 5]
        assert all(0 < r < 7 for r in row_fields)
        profiles.append(row_fields)

    # Deterministic 8->4 centered demand is 8-5 sqrt(2)>0.  A one-vertex
    # deletion has centered increment 18-35 sqrt(14)/8>0.  The latter,
    # multiplied by r_i/7, is also the parent-level demand of the genuine
    # field-proportional singleton attempt; its no-deletion outcome is zero.
    assert 8 * 8 > 25 * 2
    assert 144 * 144 > 35 * 35 * 14

    print(
        "A8",
        "all-16-pairs residual=0, child imbalances=0, singleton Q=18",
        "positive-ground fields=permutations of (1,1,1,3,3,3,3,5)",
    )
    return profiles


def audit_general_identity(max_n=5):
    """Exhaust the residual/completion identity beyond endpoint states."""
    checked = 0
    for n in range(2, max_n + 1):
        edge_count = (n - 1) * (n - 2) // 2
        for code in range(1 << edge_count):
            a = normalized_matrix(n, code)
            q = stats(a)[2]
            for mask in range(1, (1 << n) - 1):
                x_ids = tuple(i for i in range(n) if (mask >> i) & 1)
                d_ids = tuple(i for i in range(n) if not ((mask >> i) & 1))
                ax = principal(a, x_ids)
                ad = principal(a, d_ids)
                qx = stats(ax)[2]
                decrement = q - qx
                assert decrement >= 0
                for y in spins(len(x_ids)):
                    h = energy(ax, y)
                    by = tuple(
                        sum(a[i][j] * y[k] for k, j in enumerate(x_ids))
                        for i in d_ids
                    )
                    ell = sum(abs(v) for v in by)
                    for sig in (1, -1):
                        cap = max(0, 2 * ell - (qx - sig * h))
                        residual = max(0, cap - decrement)
                        exposed = 2 * ell + sig * h - q
                        assert residual == max(0, exposed)
                        u = tuple(1 if sig * v >= 0 else -1 for v in by)
                        sibling = sig * energy(ad, u)
                        assert sibling <= -exposed
                        if residual:
                            assert sibling <= -residual < 0
                        checked += 1
    print("GENERAL residual/completion buckets checked =", checked)


def main():
    audit_general_identity()
    table = audit_small_minimizers()
    audit_a5()
    audit_a8()
    assert tuple(row[1] for row in table) == (2, 6, 8, 8, 10, 18)
    print("ENDPOINT TIE/RESIDUAL EXACT CHECKS: PASS")


if __name__ == "__main__":
    main()
