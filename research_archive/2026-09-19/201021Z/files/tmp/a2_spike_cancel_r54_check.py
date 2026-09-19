#!/usr/bin/env python3
"""Exact checks for the Wave 54 A^2 spike-cancellation audit.

The order-22 construction is a complete signing with competitive-order cut
cap, but it is NOT asserted to be an exact minimizer.  A9 and A10 are the
stored exact minimizers (A10 is the stored deterministic representative).
"""

from __future__ import annotations

import itertools
import sys
from fractions import Fraction

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from box_triple_r50_check import A as A10
from envelope_block_cover_r27 import A6, A9


def projective_spins(n: int):
    for tail in itertools.product((-1, 1), repeat=n - 1):
        yield np.asarray((1,) + tail, dtype=np.int64)


def qnorm(a: np.ndarray) -> int:
    return max(abs(int(x @ a @ x)) for x in projective_spins(len(a)))


def is_signing(a: np.ndarray) -> bool:
    return (
        a.ndim == 2
        and a.shape[0] == a.shape[1]
        and np.array_equal(a, a.T)
        and np.all(np.diag(a) == 0)
        and set(map(int, a[~np.eye(len(a), dtype=bool)])) <= {-1, 1}
    )


def sylvester(order: int) -> np.ndarray:
    assert order >= 1 and order & (order - 1) == 0
    h = np.ones((1, 1), dtype=np.int64)
    while len(h) < order:
        h = np.block([[h, h], [h, -h]])
    return h


def hub_adjoin(core: np.ndarray, hubs: int) -> np.ndarray:
    r = len(core)
    d = np.ones((r + hubs, r + hubs), dtype=np.int64)
    np.fill_diagonal(d, 0)
    d[hubs:, hubs:] = core
    return d


def pair_duplicate(d: np.ndarray) -> np.ndarray:
    """Twin edge +1; every four-edge inter-pair block equals d_ij."""
    s = len(d)
    p = np.empty((2 * s, 2 * s), dtype=np.int64)
    np.fill_diagonal(p, 0)
    for i in range(s):
        p[i, s + i] = p[s + i, i] = 1
        for j in range(s):
            if i == j:
                continue
            for ii in (i, s + i):
                for jj in (j, s + j):
                    p[ii, jj] = d[i, j]
    return p


def finite_hidden_spike_construction():
    # s=8, k=6.  A6 already has positive ground 1, Q(A6)=10, and
    # A6^2=5I.  Two positive hubs create the spiky order-eight D.
    s, k = 8, 6
    assert qnorm(A6) == 10
    assert int(np.ones(6, dtype=np.int64) @ A6 @ np.ones(6, dtype=np.int64)) == 10
    assert np.array_equal(A6 @ A6, 5 * np.eye(6, dtype=np.int64))
    d0 = hub_adjoin(A6, 2)
    p = pair_duplicate(d0)
    h0 = sylvester(s)[:k]
    b = np.concatenate((h0, -h0), axis=1)
    c = A6.copy()
    a = np.block([[p, b.T], [b, c]])
    assert all(is_signing(x) for x in (d0, p, c, a))

    # Exact cap identities for the hub join and pair duplication.
    qd = qnorm(d0)
    qp = qnorm(p)
    qc = qnorm(c)
    assert qd == 10 + 2 * 2 * 6 + 2 * 1 == 36
    assert qp == 4 * qd + 2 * s == 160
    assert qc == 10

    m, n = 2 * s, 2 * s + k
    y = np.ones(m, dtype=np.int64)
    fields = p @ y
    assert np.array_equal(fields[:s], fields[s:])
    assert int(fields.max()) == 2 * s - 1
    assert np.count_nonzero(fields == 2 * s - 1) >= 4
    assert np.array_equal(b @ b.T, m * np.eye(k, dtype=np.int64))
    assert np.array_equal(b @ y, np.zeros(k, dtype=np.int64))
    assert np.array_equal(b @ fields, np.zeros(k, dtype=np.int64))

    a2 = a @ a
    cross = b @ y
    dd = a2[m:, :m] @ y
    ht = a2[m:, m:] - (n - 1) * np.eye(k, dtype=np.int64)
    internal = int(fields @ fields)
    mu = internal + int(cross @ cross) + k * (n - 1)
    rows = []
    for ww in itertools.product((-1, 1), repeat=k):
        w = np.asarray(ww, dtype=np.int64)
        x = np.concatenate((y, w))
        rows.append(int((a @ x) @ (a @ x)))
    assert np.array_equal(dd, np.zeros(k, dtype=np.int64))
    assert np.array_equal(ht, np.zeros((k, k), dtype=np.int64))
    assert set(rows) == {mu}
    assert internal == 1840 and mu == 1966

    # Analytic competitive-order bound: Q(A) <= Q(P)+Q(C)+2m sqrt(k).
    # We record its squared integer comparison rather than round the radical.
    excess_bound_sq = (2 * m) ** 2 * k
    assert excess_bound_sq == 6144
    print("scalable-template finite realization (not exact-minimal)")
    print(
        f"  n={n} m={m} k={k} Q(D)={qd} Q(P)={qp} Q(C)={qc} "
        f"Q(A)-upper={qp + qc}+sqrt({excess_bound_sq})"
    )
    print(
        f"  rmax={int(fields.max())} I={internal} X={int(cross@cross)} "
        f"||d||^2={int(dd@dd)} ||H_T||_F^2={int((ht*ht).sum())} "
        f"V=mean={mu}"
    )


def child_grounds(a: np.ndarray, m: int):
    n = len(a)
    a2 = a @ a
    for ss in itertools.combinations(range(n), m):
        s = np.asarray(ss, dtype=int)
        t = np.asarray([i for i in range(n) if i not in ss], dtype=int)
        p = a[np.ix_(s, s)]
        vals = [(y, int(y @ p @ y)) for y in projective_spins(m)]
        qs = max(abs(raw) for _, raw in vals)
        for y, raw in vals:
            if abs(raw) != qs:
                continue
            sigma = 1 if raw > 0 else -1
            r = sigma * y * (p @ y)
            b = a[np.ix_(t, s)]
            c = a[np.ix_(t, t)]
            u = b @ y
            dd = a2[np.ix_(t, s)] @ y
            ht = a2[np.ix_(t, t)] - (n - 1) * np.eye(n - m, dtype=np.int64)
            internal = int(r @ r)
            cross = int(u @ u)
            mu = internal + cross + (n - m) * (n - 1)
            rows = []
            for ww in itertools.product((-1, 1), repeat=n-m):
                w = np.asarray(ww, dtype=np.int64)
                rows.append(mu + int(2 * dd @ w + w @ ht @ w))

            # Exact spike-column identities and rational product-bias bounds.
            for loc in range(m):
                ci = sigma * int(y[loc]) * b[:, loc]
                ai = int(ci @ dd)
                gamma = sum(
                    int(y[loc] * y[j] * r[j] * (b[:, loc] @ b[:, j]))
                    for j in range(m) if j != loc
                )
                exterior = int(sigma * y[loc] * (b[:, loc] @ c @ u))
                assert ai == (n - m) * int(r[loc]) + gamma + exterior
                hi = int(ci @ ht @ ci)
                hi2 = (
                    (n-m) ** 2
                    + sum(int(b[:,loc] @ b[:,j]) ** 2 for j in range(m) if j != loc)
                    + int((c @ b[:,loc]) @ (c @ b[:,loc]))
                    - (n-m) * (n-1)
                )
                assert hi == hi2
                for theta in (Fraction(0), Fraction(1,3), Fraction(2,3), Fraction(1)):
                    expected = Fraction(mu) - 2 * theta * ai + theta * theta * hi
                    assert Fraction(min(rows)) <= expected

            yield dict(
                S=tuple(map(int, s)), y=tuple(map(int, y)), sigma=sigma,
                r=tuple(map(int, r)), rmax=int(r.max()), I=internal, X=cross,
                d=tuple(map(int, dd)), H=ht, mu=mu, V=min(rows),
            )


def finite_exact_minimizer_audit():
    a9 = list(child_grounds(A9, 8))
    hidden9 = next(z for z in a9 if z["rmax"] == 7 and not any(z["d"])
                   and not np.any(z["H"]))
    assert hidden9["I"] == 120 and hidden9["X"] == 0
    assert hidden9["V"] == hidden9["mu"] == 128

    a10m8 = list(child_grounds(A10, 8))
    hidden10 = next(z for z in a10m8 if z["rmax"] == 7 and not any(z["d"])
                    and not np.any(z["H"]))
    assert hidden10["I"] == 120 and hidden10["X"] == 0
    assert hidden10["V"] == hidden10["mu"] == 138

    a10m6 = list(child_grounds(A10, 6))
    best10 = min(a10m6, key=lambda z: z["V"])
    assert (best10["I"], best10["X"], best10["mu"], best10["V"]) == (30, 8, 74, 10)
    assert not np.any(best10["H"])
    assert sum(abs(v) for v in best10["d"]) == 32
    assert best10["V"] == best10["mu"] - 2 * sum(abs(v) for v in best10["d"])

    print("stored exact-minimizer finite audit")
    print(
        "  A9,m=8 hidden spike:",
        (hidden9["S"], hidden9["y"], hidden9["r"], hidden9["I"],
         hidden9["X"], hidden9["d"], int((hidden9["H"]**2).sum()),
         hidden9["mu"], hidden9["V"]),
    )
    print(
        "  A10,m=8 hidden spike:",
        (hidden10["S"], hidden10["y"], hidden10["r"], hidden10["I"],
         hidden10["X"], hidden10["d"], int((hidden10["H"]**2).sum()),
         hidden10["mu"], hidden10["V"]),
    )
    print(
        "  A10,m=6 best box (pure signed-d cancellation):",
        (best10["S"], best10["y"], best10["r"], best10["I"],
         best10["X"], best10["d"], int((best10["H"]**2).sum()),
         best10["mu"], best10["V"]),
    )


def main():
    finite_hidden_spike_construction()
    finite_exact_minimizer_audit()
    print("PASS exact block identities, spike certificate, and finite audits")


if __name__ == "__main__":
    main()
