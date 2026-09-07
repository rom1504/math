#!/usr/bin/env python3
"""Exact checks for Wave 47C tight-decomposition separation.

The script checks the flip-game identities, constructs an abstract block-game
obstruction, and tests common-parent and nested-chain strengthenings on exact
minimizers.  All substantive calculations are integer/exact.
"""

from __future__ import annotations

import itertools
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from tropical_decomposition_r46_check import (  # noqa: E402
    edges,
    enumerate_normalized_minimizers,
    spins,
)
from endpoint_transport_r25 import A8, A9  # noqa: E402


def edge_features(m: int) -> np.ndarray:
    """All oriented child cut features v_e=sigma*x_i*x_j."""
    x = spins(m).astype(np.int16)
    p = np.asarray([[xx[i] * xx[j] for i, j in edges(m)] for xx in x], dtype=np.int16)
    return np.concatenate((p, -p), axis=0)


def abstract_antipodal_obstruction() -> dict[str, object]:
    """A block-game profile with no tight parent/child state.

    The current block is the all-positive triangle.  Two antipodal nonground
    child features are active; all other states have positive slack.
    """
    m = 3
    v = edge_features(m)
    a = np.ones(len(edges(m)), dtype=np.int16)
    c_a = 2 * (v @ a)
    qstar = int(np.max(c_a))
    # Pick a feature of current-block energy +2 and its antipode.
    i = int(np.flatnonzero(c_a == 2)[0])
    j = int(np.flatnonzero(np.all(v == -v[i], axis=1))[0])
    active = {i, j}
    slack = np.full(len(v), 40, dtype=np.int64)
    slack[list(active)] = 0
    # Keep the abstract completion profile nonnegative, as every genuine
    # outside polynomial has maximum at least its (zero) uniform mean.
    q = 100
    k = q - c_a - slack

    rows = []
    gamma_bool = -10**9
    for fbits in itertools.product((0, 1), repeat=len(a)):
        f = np.asarray(fbits, dtype=np.int16)
        margins = slack + 4 * ((v * a[None, :]) @ f)
        margin = int(np.min(margins))
        gamma_bool = max(gamma_bool, margin)
        b = a * (1 - 2 * f)
        new_payoff = k + 2 * (v @ b)
        assert int(np.max(new_payoff)) >= q
        assert margin == q - int(np.max(new_payoff))
        rows.append((tuple(int(x) for x in f), margin))
    assert gamma_bool == 0
    assert all(c_a[t] < qstar for t in active)

    # The half-half active distribution is the exact fractional-dual
    # certificate: E s=0 and every E[a_e v_e]=0.
    mu = (v[i] + v[j]) / 2
    assert np.all(mu == 0)
    dual_value = 0.0
    return {
        "m": m,
        "a": tuple(int(x) for x in a),
        "active_features": [tuple(int(x) for x in v[t]) for t in (i, j)],
        "active_child_energies": [int(c_a[t]) for t in (i, j)],
        "q_star": qstar,
        "gamma_bool": gamma_bool,
        "fractional_dual_value": dual_value,
        "flip_margins": rows,
    }


def abstract_profile_lowering_obstruction() -> dict[str, object]:
    """Failure plus block minimality need not permit a profile-lowering move.

    This five-vertex block has norm 12 although the order-five optimum is 8.
    An antipodal active pair makes a cap-preserving flip possible only when
    its signed flip sum is zero.  Exact exhaustion shows that no such safe
    flip lowers the block norm.
    """
    m = 5
    v = edge_features(m).astype(np.int64)
    a = np.asarray((-1, -1, -1, -1, -1, -1, -1, -1, 1, 1), dtype=np.int64)
    ca = 2 * (v @ a)
    qstar = int(np.max(ca))
    assert qstar == 12
    target = np.asarray((-1, -1, -1, -1, 1, 1, 1, 1, 1, 1), dtype=np.int64)
    i = int(np.flatnonzero(np.all(v == target, axis=1))[0])
    j = int(np.flatnonzero(np.all(v == -target, axis=1))[0])
    assert int(ca[i]) == 4 and int(ca[j]) == -4
    slack = np.full(len(v), 100, dtype=np.int64)
    slack[[i, j]] = 0

    min_child_norm = 10**9
    safe_flips = []
    safe_lower = []
    game_margin = -10**9
    for fbits in itertools.product((0, 1), repeat=len(a)):
        f = np.asarray(fbits, dtype=np.int64)
        b = a * (1 - 2 * f)
        qb = int(np.max(2 * (v @ b)))
        min_child_norm = min(min_child_norm, qb)
        margins = slack + 4 * ((v * a[None, :]) @ f)
        margin = int(np.min(margins))
        game_margin = max(game_margin, margin)
        if margin >= 0:
            safe_flips.append((tuple(int(z) for z in f), qb))
            if qb < qstar:
                safe_lower.append((tuple(int(z) for z in f), qb))
    assert min_child_norm == 8
    assert game_margin == 0
    assert safe_flips and not safe_lower
    assert all(qb >= qstar for _, qb in safe_flips)
    return {
        "m": m,
        "a": tuple(int(z) for z in a),
        "q_star": qstar,
        "minimum_child_norm": min_child_norm,
        "active_features": [tuple(int(z) for z in v[t]) for t in (i, j)],
        "active_child_energies": [int(ca[t]) for t in (i, j)],
        "cap_preserving_flip_count": len(safe_flips),
        "cap_preserving_lowering_count": len(safe_lower),
        "cap_preserving_child_norms": sorted(set(qb for _, qb in safe_flips)),
    }


def subset_data(A: np.ndarray) -> dict[str, object]:
    """Exact maximal-principal and inherited-ground data for one signing."""
    n = len(A)
    ee = edges(n)
    avec = np.asarray([A[i, j] for i, j in ee], dtype=np.int16)
    x = spins(n).astype(np.int16)
    prod = np.asarray([[xx[i] * xx[j] for i, j in ee] for xx in x], dtype=np.int16)
    raw = 2 * (prod @ avec)
    qn = int(np.max(np.abs(raw)))
    parent_ids = np.flatnonzero(np.abs(raw) == qn)
    sigmas = np.where(raw[parent_ids] > 0, 1, -1).astype(np.int16)

    masks = np.arange(1 << n, dtype=np.int64)
    sizes = np.asarray([bin(int(s)).count("1") for s in masks], dtype=np.int8)
    edge_masks = np.asarray(
        [[int(((s >> i) & 1) and ((s >> j) & 1)) for i, j in ee] for s in masks],
        dtype=np.int16,
    )
    internal = 2 * (prod @ (edge_masks * avec[None, :]).T)
    qmask = np.max(np.abs(internal), axis=0)
    qstar = np.asarray([int(np.max(qmask[sizes == m])) for m in range(n + 1)], dtype=np.int16)

    reports = []
    for pid, sigma in zip(parent_ids, sigmas):
        good = (qmask == qstar[sizes]) & (sigma * internal[pid] == qstar[sizes])
        by_size = {m: [int(s) for s in masks[(sizes == m) & good]] for m in range(n + 1)}
        all_m = all(by_size[m] for m in range(n + 1))

        reachable = np.zeros(1 << n, dtype=bool)
        parent = np.full(1 << n, -1, dtype=np.int64)
        reachable[0] = bool(good[0])
        first_broken = None
        for m in range(1, n + 1):
            for s in by_size[m]:
                bits = [i for i in range(n) if (s >> i) & 1]
                for i in bits:
                    p = s ^ (1 << i)
                    if reachable[p]:
                        reachable[s] = True
                        parent[s] = p
                        break
            if first_broken is None and not np.any(reachable[sizes == m]):
                first_broken = m
        nested = bool(reachable[(1 << n) - 1])

        # Accessibility and exchange of the full good-set family.
        accessibility_bad = []
        for s in np.flatnonzero(good):
            if sizes[s] <= 1:
                continue
            if not any(good[s ^ (1 << i)] for i in range(n) if (s >> i) & 1):
                accessibility_bad.append(int(s))
        exchange_bad = []
        good_ids = [int(s) for s in np.flatnonzero(good)]
        for u in good_ids:
            for w in good_ids:
                if sizes[u] <= sizes[w]:
                    continue
                candidates = [i for i in range(n) if ((u >> i) & 1) and not ((w >> i) & 1)]
                if not any(good[w | (1 << i)] for i in candidates):
                    exchange_bad.append((u, w))
                    break
            if exchange_bad:
                break

        reports.append(
            {
                "parent_spin": tuple(int(z) for z in x[pid]),
                "sigma": int(sigma),
                "all_m": all_m,
                "nested": nested,
                "first_broken_size": first_broken,
                "good_by_size": by_size,
                "reachable_by_size": {
                    m: [int(s) for s in masks[(sizes == m) & reachable]] for m in range(n + 1)
                },
                "accessibility_bad": accessibility_bad,
                "exchange_bad": exchange_bad,
            }
        )
    return {
        "n": n,
        "q_n": qn,
        "q_star": [int(z) for z in qstar],
        "parent_count": len(reports),
        "reports": reports,
    }


def summarize_one(A: np.ndarray, name: str, detailed: bool = False) -> dict[str, object]:
    z = subset_data(A)
    rr = z["reports"]
    out: dict[str, object] = {
        "name": name,
        "n": z["n"],
        "q_n": z["q_n"],
        "q_star": z["q_star"],
        "parent_count": z["parent_count"],
        "all_m_parent_count": sum(r["all_m"] for r in rr),
        "nested_parent_count": sum(r["nested"] for r in rr),
        "first_broken_histogram": dict(Counter(r["first_broken_size"] for r in rr)),
        "accessibility_bad_parent_count": sum(bool(r["accessibility_bad"]) for r in rr),
        "exchange_bad_parent_count": sum(bool(r["exchange_bad"]) for r in rr),
    }
    if detailed:
        out["reports"] = rr
    return out


def exhaust_minimizers(max_n: int = 8) -> list[dict[str, object]]:
    """Exhaust common-parent and nested-chain properties through max_n."""
    ans = []
    for n in range(2, max_n + 1):
        qn, minimizers = enumerate_normalized_minimizers(n)
        all_m_fail = 0
        nested_fail = 0
        every_parent_all_m_fail = 0
        every_parent_nested_fail = 0
        no_fixed_parent_greedoid = 0
        first_accessibility_failure = None
        first_exchange_failure = None
        first_nested_failure = None
        for ai, A in enumerate(minimizers):
            z = summarize_one(A, f"n{n}_class{ai}", detailed=True)
            assert z["q_n"] == qn
            if z["all_m_parent_count"] == 0:
                all_m_fail += 1
            if z["nested_parent_count"] == 0:
                nested_fail += 1
                if first_nested_failure is None:
                    first_nested_failure = ai
            if z["all_m_parent_count"] != z["parent_count"]:
                every_parent_all_m_fail += 1
            if z["nested_parent_count"] != z["parent_count"]:
                every_parent_nested_fail += 1
            reports = z["reports"]
            if not any(not r["accessibility_bad"] and not r["exchange_bad"] for r in reports):
                no_fixed_parent_greedoid += 1
            if first_accessibility_failure is None and any(r["accessibility_bad"] for r in reports):
                first_accessibility_failure = ai
            if first_exchange_failure is None and any(r["exchange_bad"] for r in reports):
                first_exchange_failure = ai
        ans.append(
            {
                "n": n,
                "q_n": qn,
                "classes": len(minimizers),
                "no_all_m_common_parent": all_m_fail,
                "no_nested_chain": nested_fail,
                "not_every_parent_all_m": every_parent_all_m_fail,
                "not_every_parent_nested": every_parent_nested_fail,
                "first_nested_failure_class": first_nested_failure,
                "no_fixed_parent_greedoid": no_fixed_parent_greedoid,
                "first_accessibility_failure_class": first_accessibility_failure,
                "first_exchange_failure_class": first_exchange_failure,
            }
        )
    return ans


def mask_vertices(mask: int, n: int) -> tuple[int, ...]:
    return tuple(i for i in range(n) if (mask >> i) & 1)


def main() -> None:
    obstruction = abstract_antipodal_obstruction()
    print("abstract antipodal obstruction", obstruction)
    lowering = abstract_profile_lowering_obstruction()
    print("abstract profile-lowering obstruction", lowering)

    a8 = summarize_one(A8, "A8", detailed=True)
    a9 = summarize_one(A9, "A9", detailed=True)
    print("stored summary A8", {k: v for k, v in a8.items() if k != "reports"})
    print("stored summary A9", {k: v for k, v in a9.items() if k != "reports"})

    # Give a human-readable exact wall for the first A8 parent ground.
    r = a8["reports"][0]
    broken = r["first_broken_size"]
    explicit = {
        "parent_spin": r["parent_spin"],
        "sigma": r["sigma"],
        "first_broken_size": broken,
        "reachable_previous": [mask_vertices(s, 8) for s in r["reachable_by_size"][broken - 1]],
        "all_good_current": [mask_vertices(s, 8) for s in r["good_by_size"][broken]],
        "accessibility_bad": [mask_vertices(s, 8) for s in r["accessibility_bad"][:5]],
        "exchange_bad": [
            (mask_vertices(u, 8), mask_vertices(w, 8)) for u, w in r["exchange_bad"][:5]
        ],
    }
    print("explicit A8 nested wall", explicit)

    exhaustion = exhaust_minimizers(8)
    print("exact-minimizer chain exhaustion", exhaustion)
    assert not any(z["no_all_m_common_parent"] for z in exhaustion)
    print("PASS tight_separation_r47_check")


if __name__ == "__main__":
    main()
