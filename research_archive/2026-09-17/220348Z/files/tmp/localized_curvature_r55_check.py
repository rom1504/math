#!/usr/bin/env python3
"""Checks for the Wave 55 localized-curvature replacement theorem.

The order-22 signing is the Wave 54 generic hidden-spike construction and is
not an exact minimizer.  A9 and A10 are stored exact minimizers.
"""

from __future__ import annotations

import itertools
import math
import sys

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from a2_spike_cancel_r54_check import (
    child_grounds,
    hub_adjoin,
    pair_duplicate,
    projective_spins,
    sylvester,
)
from box_triple_r50_check import A as A10
from envelope_block_cover_r27 import A6, A8, A9


def qnorm_batched(a: np.ndarray, batch: int = 1 << 17):
    """Exact projective enumeration, returning cap, raw value, and word."""
    n = len(a)
    best = -1
    best_raw = 0
    best_x = None
    for start in range(0, 1 << (n - 1), batch):
        ids = np.arange(
            start, min(start + batch, 1 << (n - 1)), dtype=np.uint32
        )
        tails = 1 - 2 * (
            (ids[:, None] >> np.arange(n - 1, dtype=np.uint32)) & 1
        ).astype(np.int16)
        x = np.concatenate(
            (np.ones((len(ids), 1), dtype=np.int16), tails), axis=1
        )
        vals = np.einsum("bi,ij,bj->b", x, a, x, optimize=True)
        loc = int(np.argmax(np.abs(vals)))
        val = abs(int(vals[loc]))
        if val > best:
            best = val
            best_raw = int(vals[loc])
            best_x = x[loc].copy()
    return best, best_raw, best_x


def top_prefix(r: np.ndarray):
    """Prefix maximizing R_h/sqrt(h), with a deterministic tie rule."""
    order = sorted(range(len(r)), key=lambda i: (-int(r[i]), i))
    sums = np.cumsum([int(r[i]) for i in order])
    scores = [int(sums[h - 1]) / math.sqrt(h) for h in range(1, len(r) + 1)]
    h = max(range(1, len(r) + 1), key=lambda j: (scores[j - 1], -j))
    return order[:h], int(sums[h - 1]), scores[h - 1]


def excess_edge_set(a: np.ndarray, rec: dict, local_u: list[int]):
    """Choose r_i child-positive incident edges for every i in U."""
    s = list(rec["S"])
    y = np.asarray(rec["y"], dtype=np.int64)
    sigma = int(rec["sigma"])
    r = np.asarray(rec["r"], dtype=np.int64)
    p = a[np.ix_(s, s)]
    edge_set = set()
    for i in local_u:
        positive = [
            j
            for j in range(len(s))
            if j != i and sigma * int(y[i]) * int(p[i, j]) * int(y[j]) == 1
        ]
        assert len(positive) == (len(s) - 1 + int(r[i])) // 2
        assert int(r[i]) <= len(positive)
        for j in positive[: int(r[i])]:
            edge_set.add(tuple(sorted((s[i], s[j]))))
    return sorted(edge_set)


def check_lorentz_and_edges(a: np.ndarray, rec: dict):
    r = np.asarray(rec["r"], dtype=np.int64)
    u, rsum, kappa = top_prefix(r)
    harmonic = sum(1.0 / j for j in range(1, len(r) + 1))
    assert int(r @ r) <= kappa * kappa * harmonic + 1e-9
    h_edges = excess_edge_set(a, rec, u)
    assert len(h_edges) >= math.ceil(rsum / 2)
    incident = set(itertools.chain.from_iterable(h_edges))
    eta = math.sqrt(8 * len(h_edges) * len(incident) * math.log(2))
    return u, rsum, h_edges, eta


def escape_profile(a: np.ndarray, h_edges: list[tuple[int, int]]):
    q, _, _ = qnorm_batched(a)
    e = len(h_edges)
    best_escape = 10**9
    ground_neg = []
    for x in projective_spins(len(a)):
        raw = int(x @ a @ x)
        for sigma in (1, -1):
            deficit = q - sigma * raw
            negative = sum(
                sigma * int(a[i, j]) * int(x[i]) * int(x[j]) < 0
                for i, j in h_edges
            )
            if negative > e / 4:
                best_escape = min(best_escape, deficit)
            if deficit == 0:
                ground_neg.append(negative)
    return q, best_escape, min(ground_neg), max(ground_neg)


def hidden_record(a: np.ndarray):
    return next(
        rec
        for rec in child_grounds(a, 8)
        if rec["rmax"] == 7 and not any(rec["d"]) and not np.any(rec["H"])
    )


def finite_exact_audit():
    for name, a in (("A9", A9), ("A10", A10)):
        rec = hidden_record(a)
        u, rsum, h_edges, eta = check_lorentz_and_edges(a, rec)
        q, best_escape, min_ground, max_ground = escape_profile(a, h_edges)
        assert rec["I"] == 120 and rec["X"] == 0
        assert best_escape == 0
        assert max_ground > len(h_edges) / 4
        print(
            f"{name} exact hidden spike: q={q} U={u} R_U={rsum} "
            f"e={len(h_edges)} eta={eta:.12f} "
            f"min-escape-deficit={best_escape} "
            f"ground-negative-range=[{min_ground},{max_ground}]"
        )


def finite_generic_replacement():
    # Wave 54 order-22 hidden-spike template.
    s, k = 8, 6
    w0 = sylvester(s)[:k]
    b = np.concatenate((w0, -w0), axis=1)
    c = A6.copy()

    d_hub = hub_adjoin(A6, 2)
    p_hub = pair_duplicate(d_hub)
    a_hub = np.block([[p_hub, b.T], [b, c]])

    # Concrete internal-S replacement: replace the hub-adjoined D by the
    # stored exact A8 before pair duplication; B and C remain unchanged.
    p_repl = pair_duplicate(A8)
    a_repl = np.block([[p_repl, b.T], [b, c]])

    q_hub, raw_hub, x_hub = qnorm_batched(a_hub)
    q_repl, raw_repl, x_repl = qnorm_batched(a_repl)
    assert (q_hub, q_repl) == (170, 118)
    assert q_repl < q_hub

    y = np.ones(2 * s, dtype=np.int64)
    fields = p_hub @ y
    assert int(fields @ fields) == 1840
    assert np.array_equal(b @ y, np.zeros(k, dtype=np.int64))
    assert np.array_equal(b @ fields, np.zeros(k, dtype=np.int64))

    print(
        "order-22 generic template (not exact-minimal): "
        f"I=1840 X=0, internal replacement cap {q_hub}->{q_repl}"
    )
    print(
        "  old maximizer=",
        raw_hub,
        "".join("+" if z == 1 else "-" for z in x_hub),
    )
    print(
        "  new maximizer=",
        raw_repl,
        "".join("+" if z == 1 else "-" for z in x_repl),
    )


def attenuation_constant_audit():
    # Pure scalar audit of the two response branches.  With
    # z=(1-t/(4e))a_H, the fractional energy loss is t*S/(2e).
    for e in (17, 100):
        for lam in (0.1, 0.5, 1.0):
            for t in (0.25 * e, 2.0 * e, 4.0 * e):
                low_gain = lam * t / 2
                high_gain = t - t / 2
                assert min(low_gain, high_gain) == lam * t / 2
                coefficient = 1 - t / (4 * e)
                assert 0 <= coefficient <= 1
    print("attenuation constants: PASS (domain 0<t<=4e)")


def main():
    attenuation_constant_audit()
    finite_exact_audit()
    finite_generic_replacement()
    print("PASS localized replacement, Lorentz-prefix, and finite diagnostics")


if __name__ == "__main__":
    main()
