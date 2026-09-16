#!/usr/bin/env python3
"""Exact finite audit for the Wave 58 two-column/switching memo.

The only floating quantity is (m/n)^(3/2) in the favorable threshold.  Every
decision has a reported positive margin; counts, rows, deficits, LP values,
and supporting-line certificates are otherwise checked with integers or
Fractions.
"""

from __future__ import annotations

import itertools
import math
import sys
from collections import Counter
from fractions import Fraction

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from envelope_block_cover_r27 import A8, A9, projective_spins


def qnorm(a: np.ndarray) -> int:
    return max(abs(int(x @ a @ x)) for x in projective_spins(len(a)))


def normalize(x: np.ndarray) -> tuple[int, ...]:
    if x[0] < 0:
        x = -x
    return tuple(int(z) for z in x)


def build_columns(a: np.ndarray, m: int, t: float = 0.0):
    n = len(a)
    q = qnorm(a)
    selectors = list(itertools.combinations(range(n), m))
    selector_arrays = [np.asarray(s, dtype=int) for s in selectors]
    child_q = [qnorm(a[np.ix_(s, s)]) for s in selectors]
    p2 = Fraction(m * (m - 1), n * (n - 1))
    p2f = float(p2)
    b_allow = ((m / n) ** 1.5 - p2f) * q
    records = {}
    decision_margin = math.inf
    for x in projective_spins(n):
        key = tuple(int(z) for z in x)
        raw = int(x @ a @ x)
        row = int((a @ x) @ (a @ x))
        for sigma in (-1, 1):
            energy = sigma * raw
            deficit = q - energy
            gs = []
            local_energies = []
            for s, idx, qs in zip(selectors, selector_arrays, child_q):
                local = sigma * int(x[idx] @ a[np.ix_(s, s)] @ x[idx])
                g = qs - local - p2f * deficit - b_allow
                decision_margin = min(decision_margin, abs(g - t))
                gs.append(g)
                local_energies.append(local)
                # Retained-deficit decomposition (10.793), pointwise.
                centered = local - p2f * energy
                local_deficit = qs - local
                lhs = qs - p2f * q
                rhs = centered + local_deficit - p2f * deficit
                assert abs(lhs - rhs) < 1e-11
            count = sum(g <= t for g in gs)
            records[(key, sigma)] = {
                "x": x.copy(),
                "sigma": sigma,
                "energy": energy,
                "delta": deficit,
                "row": row,
                "count": count,
                "g": np.asarray(gs),
                "local": np.asarray(local_energies, dtype=int),
            }
    assert decision_margin > 1e-9
    return {
        "a": a,
        "n": n,
        "m": m,
        "q": q,
        "N": len(selectors),
        "selectors": selectors,
        "child_q": child_q,
        "p2": p2,
        "b": b_allow,
        "records": records,
        "margin": decision_margin,
    }


def orientation_and_switching_identities(data) -> None:
    a = data["a"]
    n = data["n"]
    q = data["q"]
    p = data["m"] / n
    p2 = float(data["p2"])
    records = data["records"]
    child_q = np.asarray(data["child_q"], dtype=float)

    # Orientation sum (R58.10), on every oriented cut and selector.
    target = 2.0 * (child_q - p ** 1.5 * q)
    for (key, sigma), rec in records.items():
        mate = records[(key, -sigma)]
        assert np.max(np.abs(rec["g"] + mate["g"] - target)) < 2e-11

    # Singleton switching and row identities on every state.
    a2 = a @ a
    for (key, sigma), rec in records.items():
        x = np.asarray(key, dtype=np.int64)
        row_sum = 0
        for i in range(n):
            y = x.copy()
            y[i] *= -1
            switched = records[(normalize(y), sigma)]
            w = np.asarray([
                sum(
                    int(a[i, j] * x[i] * x[j])
                    for j in s if j != i
                ) if i in s else 0
                for s in data["selectors"]
            ], dtype=float)
            parent_w = sum(int(a[i, j] * x[i] * x[j])
                           for j in range(n) if j != i)
            predicted = 4 * sigma * (w - p2 * parent_w)
            assert np.max(np.abs(switched["g"] - rec["g"] - predicted)) < 2e-11
            row_increment = switched["row"] - rec["row"]
            predicted_row = 4 * ((n - 1) - int(x[i] * (a2 @ x)[i]))
            assert row_increment == predicted_row
            row_sum += row_increment
        assert row_sum == 4 * (n * (n - 1) - rec["row"])


def objective(rec, n_selectors: int, price: Fraction, logarithmic: bool):
    if rec["count"] == 0:
        return math.inf
    if logarithmic:
        return math.log(n_selectors / rec["count"]) + float(price * rec["row"])
    return n_selectors / rec["count"] + float(price * rec["row"])


def a8_full_bare_barrier(data) -> dict:
    assert (data["q"], data["N"]) == (20, 56)
    records = data["records"]
    key = ((1, -1, -1, -1, -1, -1, 1, 1), -1)
    rec = records[key]
    assert (rec["count"], rec["row"], rec["delta"]) == (32, 72, 20)
    mate = records[(key[0], 1)]
    assert (mate["count"], mate["row"], mate["delta"]) == (32, 72, 20)

    x = np.asarray(key[0], dtype=np.int64)
    neighbors = []
    for i in range(8):
        y = x.copy()
        y[i] *= -1
        neighbors.append(records[(normalize(y), -1)])
    assert Counter((z["count"], z["row"], z["delta"]) for z in neighbors) == Counter({
        (19, 64, 8): 4,
        (15, 64, 32): 4,
    })

    price = Fraction(1, 20)
    for logarithmic in (True, False):
        base = objective(rec, 56, price, logarithmic)
        assert objective(mate, 56, price, logarithmic) == base
        assert all(objective(z, 56, price, logarithmic) > base for z in neighbors)
        global_value = min(objective(z, 56, price, logarithmic)
                           for z in records.values())
        global_types = {
            (z["count"], z["row"], z["delta"])
            for z in records.values()
            if abs(objective(z, 56, price, logarithmic) - global_value) < 1e-12
        }
        assert global_types == {(36, 8, 20)}
        assert global_value < base

    # Exact positivity of the best-neighbor increments in (R58.17).
    assert math.log(Fraction(32, 19)) > Fraction(2, 5)
    assert Fraction(56, 19) - Fraction(7, 4) - Fraction(2, 5) > 0

    # Best block switches by projective distance.  Both scalarizations choose
    # the displayed type at each distance.
    expected = [
        (32, 72, 20),
        (19, 64, 8),
        (20, 40, 4),
        (39, 32, 16),
        (36, 8, 20),
    ]
    block_table = []
    for k in range(5):
        candidates = []
        for block in itertools.combinations(range(8), k):
            y = x.copy()
            y[list(block)] *= -1
            z = records[(normalize(y), -1)]
            if z["count"]:
                candidates.append((block, z))
        best_log = min(candidates, key=lambda bz: objective(bz[1], 56, price, True))
        best_exp = min(candidates, key=lambda bz: objective(bz[1], 56, price, False))
        typ_log = (best_log[1]["count"], best_log[1]["row"], best_log[1]["delta"])
        typ_exp = (best_exp[1]["count"], best_exp[1]["row"], best_exp[1]["delta"])
        assert typ_log == typ_exp == expected[k]
        block_table.append((k, typ_log, best_log[0]))
    return {"neighbor_types": Counter((z["count"], z["row"], z["delta"])
                                      for z in neighbors),
            "block_table": block_table}


def finite_types(data, sector_cap=None):
    ans = set()
    for rec in data["records"].values():
        if rec["count"] == 0:
            continue
        if sector_cap is not None and rec["delta"] > sector_cap:
            continue
        ans.add((rec["count"], rec["row"], rec["delta"]))
    return sorted(ans)


def exact_row_lp(types, n_selectors: int, row_cap: int):
    """Enumerate all one- and two-point captured laws exactly."""
    best = None
    supports = []
    for typ in types:
        count, row, _ = typ
        if row <= row_cap:
            value = Fraction(n_selectors, count)
            if best is None or value < best:
                best, supports = value, [((Fraction(1), typ),)]
            elif value == best:
                supports.append(((Fraction(1), typ),))
    for low in types:
        for high in types:
            if not (low[1] < row_cap < high[1]):
                continue
            w_low = Fraction(high[1] - row_cap, high[1] - low[1])
            w_high = 1 - w_low
            value = (w_low * Fraction(n_selectors, low[0])
                     + w_high * Fraction(n_selectors, high[0]))
            if best is None or value < best:
                best = value
                supports = [((w_low, low), (w_high, high))]
            elif value == best:
                supports.append(((w_low, low), (w_high, high)))
    assert best is not None
    return best, supports


def strict_cap_checks(a8, a9) -> dict:
    types8 = finite_types(a8)
    m8, supports8 = exact_row_lp(types8, 56, 40)
    assert m8 == Fraction(371, 260)
    assert any(
        set((w, typ[:2]) for w, typ in support)
        == {(Fraction(3, 4), (39, 32)), (Fraction(1, 4), (40, 64))}
        for support in supports8
    )
    assert max(count for count, row, _ in types8 if row <= 40) == 39

    # Exact supporting-line certificate in the (R,e^h) plane.
    gamma8 = Fraction(7, 6240)
    support_value8 = Fraction(287, 195)
    assert min(Fraction(56, c) + gamma8 * r for c, r, _ in types8) == support_value8
    active8 = {(c, r, d) for c, r, d in types8
               if Fraction(56, c) + gamma8 * r == support_value8}
    assert active8 == {(39, 32, 16), (40, 64, 40)}
    assert support_value8 - gamma8 * 40 == m8
    # Captured-to-reference conversion (R58.1).
    assert Fraction(260, 371) * Fraction(1, 4) * Fraction(7, 5) == Fraction(91, 371)
    assert Fraction(260, 371) * Fraction(3, 4) * Fraction(56, 39) == Fraction(280, 371)
    assert Fraction(260, 371) > Fraction(39, 56)

    # Correct-sector restriction removes the A8 strict-cap improvement.
    correct8 = finite_types(a8, sector_cap=a8["q"])
    m8_correct, _ = exact_row_lp(correct8, 56, 40)
    assert m8_correct == Fraction(56, 39)

    # A9 has a genuine correct-sector, strict-cap, adjacent switching pair.
    correct9 = finite_types(a9, sector_cap=a9["q"])
    m9, supports9 = exact_row_lp(correct9, 84, 64)
    assert m9 == Fraction(217, 80)
    target9 = {(Fraction(1, 2), (30, 56, 24)),
               (Fraction(1, 2), (32, 72, 24))}
    assert any(set(support) == target9 for support in supports9)
    gamma9 = Fraction(7, 640)
    support_value9 = Fraction(273, 80)
    assert min(Fraction(84, c) + gamma9 * r for c, r, _ in correct9) == support_value9
    assert support_value9 - gamma9 * 64 == m9
    assert max(c for c, r, _ in correct9 if r <= 64) == 30

    low_keys = [key for key, rec in a9["records"].items()
                if (rec["count"], rec["row"], rec["delta"]) == (30, 56, 24)]
    high_keys = [key for key, rec in a9["records"].items()
                 if (rec["count"], rec["row"], rec["delta"]) == (32, 72, 24)]
    adjacent = False
    for (x, sigma) in low_keys:
        xx = np.asarray(x)
        for (y, tau) in high_keys:
            if sigma != tau:
                continue
            yy = np.asarray(y)
            distance = min(int(np.sum(xx != yy)), int(np.sum(xx != -yy)))
            adjacent |= distance == 1
    assert adjacent

    # Jensen is strict for the A8 captured pair.
    eh = 0.25 * math.log(7 / 5) + 0.75 * math.log(56 / 39)
    selected_surprise = math.log(371 / 260)
    assert selected_surprise > eh
    return {
        "A8_full": {"M": str(m8), "Z": str(1 / m8),
                    "best_strict_single": "39/56"},
        "A8_correct_sector": {"M": str(m8_correct)},
        "A9_correct_sector": {"M": str(m9), "Z": str(1 / m9),
                              "adjacent_pair": True},
        "Jensen_gap_A8": selected_surprise - eh,
    }


def constant_slack_check(data) -> None:
    """Finite exact audit of (R58.7) with b=2 across row caps."""
    types = finite_types(data)
    for row_cap in sorted({r for _, r, _ in types if r > 0}):
        m_value, _ = exact_row_lp(types, data["N"], row_cap)
        z_value = 1 / m_value
        relaxed = max(Fraction(c, data["N"]) for c, r, _ in types
                      if r <= 2 * row_cap)
        assert z_value <= 2 * relaxed


def main() -> None:
    a8 = build_columns(A8, 5)
    a9 = build_columns(A9, 6)
    orientation_and_switching_identities(a8)
    barrier = a8_full_bare_barrier(a8)
    strict = strict_cap_checks(a8, a9)
    constant_slack_check(a8)
    constant_slack_check(a9)
    print({
        "decision_margins": {"A8": a8["margin"], "A9": a9["margin"]},
        "A8_barrier": barrier,
        "strict_cap": strict,
    })
    print("PASS two-column KKT, retained-deficit switching, and full-bare local barrier")


if __name__ == "__main__":
    main()
