#!/usr/bin/env python3
"""Exact checks for Wave 36 simultaneous replacement attack.

All paths and generated data remain under the repository-local tmp directory.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product
from math import floor
from typing import Optional

import numpy as np
from scipy.optimize import linprog

from coherent_signatures_r29_search import A6, A8, A9, x_of


def states(n: int):
    for sigma in (-1, 1):
        for tail in product((-1, 1), repeat=n - 1):
            yield sigma, np.asarray((1,) + tail, dtype=np.int64)


def qnorm(a: np.ndarray) -> int:
    return max(int(sigma * x @ a @ x) for sigma, x in states(len(a)))


def response_table(a: np.ndarray, selector: tuple[int, ...]):
    """Return canonical projective labels, child energies, and outside widths."""
    n = len(a)
    s = np.asarray(selector, dtype=np.int64)
    outside = np.asarray([i for i in range(n) if i not in selector], dtype=np.int64)
    child = a[np.ix_(s, s)]
    cross = a[np.ix_(s, outside)]
    outer = a[np.ix_(outside, outside)]
    labels = [np.asarray((1,) + tail, dtype=np.int64)
              for tail in product((-1, 1), repeat=len(selector) - 1)]
    rows = []
    for y in labels:
        energy = int(y @ child @ y)
        z = {}
        for sigma in (-1, 1):
            z[sigma] = max(
                int(sigma * (x @ outer @ x + 2 * y @ cross @ x))
                for x in map(np.asarray, product((-1, 1), repeat=len(outside)))
            )
        assert (z[-1] + z[1]) % 2 == 0
        rows.append({
            "label": y,
            "energy": energy,
            "width": (z[-1] + z[1]) // 2,
            "zminus": z[-1],
            "zplus": z[1],
        })
    return rows


def favorable_rows(a: np.ndarray, selector: tuple[int, ...], allowance: int):
    rows = response_table(a, selector)
    child_q = max(abs(row["energy"]) for row in rows)
    return [row for row in rows if child_q - abs(row["energy"]) <= allowance]


def hard_favorable_value(a: np.ndarray, selector: tuple[int, ...], allowance: int):
    """Exact continuous one-block projective-doubled favorable value.

    A paired pair of response orientations supplies the exact lower bound
    ``max_y (Z_+(y)+Z_-(y))/2``.  The LP is used only to find a small rational
    fractional block certifying the reverse inequality; every resulting
    inequality is then checked with ``Fraction`` arithmetic.
    """
    rows = favorable_rows(a, selector, allowance)
    edges = list(combinations(range(len(selector)), 2))
    drows, zrows = [], []
    for row in rows:
        y = row["label"]
        base = np.asarray([y[i] * y[j] for i, j in edges], dtype=float)
        for sigma in (-1, 1):
            drows.append(sigma * base)
            zrows.append(row["zminus"] if sigma == -1 else row["zplus"])
    dmat = np.asarray(drows, dtype=float)
    zvec = np.asarray(zrows, dtype=float)
    # min_{xi,t} t subject to Z_u+2<xi,d_u> <= t.
    lp = linprog(
        np.r_[np.zeros(len(edges)), 1.0],
        A_ub=np.c_[2.0 * dmat, -np.ones(len(zvec))],
        b_ub=-zvec,
        bounds=[(-1.0, 1.0)] * len(edges) + [(None, None)],
        method="highs",
    )
    assert lp.success
    lower = max((row["zminus"] + row["zplus"]) // 2 for row in rows)
    assert abs(lp.fun - lower) < 1e-8
    xi = [Fraction(float(value)).limit_denominator(10_000)
          for value in lp.x[:-1]]
    assert all(-1 <= value <= 1 for value in xi)
    for drow, zvalue in zip(drows, zrows):
        score = Fraction(int(zvalue)) + 2 * sum(
            value * int(direction) for value, direction in zip(xi, drow)
        )
        assert score <= lower
    return lower


def audit_widths(a: np.ndarray, name: str, m: int):
    n = len(a)
    q = qnorm(a)
    allowance = floor(((m / n) ** 1.5 - (m / n) ** 2) * q + 1e-12)
    selectors = list(combinations(range(n), m))
    full_rows = [int((a @ x) @ (a @ x))
                 for x in (x_of(mask, n) for mask in range(1 << (n - 1)))]
    assert max(full_rows) <= 2 * n * (n - 1)
    data = []
    for s in selectors:
        favorable = favorable_rows(a, s, allowance)
        assert favorable
        widths = [row["width"] for row in favorable]
        # Parent optimality alone gives width <= q for every local label.
        assert max(widths) <= q
        hard_value = hard_favorable_value(a, s, allowance)
        data.append((s, q - max(widths), len(favorable), q - hard_value))
    return {
        "name": name,
        "n": n,
        "m": m,
        "q": q,
        "allowance": allowance,
        "width_deficits": sorted({v: sum(d == v for _, d, _, _ in data)
                                   for v in {d for _, d, _, _ in data}}.items()),
        "hard_game_deficits": sorted({v: sum(g == v for _, _, _, g in data)
                                       for v in {g for _, _, _, g in data}}.items()),
        "favorable_count_range": (min(c for _, _, c, _ in data),
                                   max(c for _, _, c, _ in data)),
        "mean_best_width_deficit": Fraction(sum(d for _, d, _, _ in data), len(data)),
        "mean_hard_game_deficit": Fraction(sum(g for _, _, _, g in data), len(data)),
        "full_row_range": (min(full_rows), max(full_rows)),
    }


def anchor_orient(label: np.ndarray, selector: tuple[int, ...], anchor: int):
    value = int(label[selector.index(anchor)])
    return label * value


def corrected_conflict(selectors, labels, anchor: int):
    count = len(selectors)
    out = Fraction(0)
    for i in range(max(max(s) for s in selectors) + 1):
        if i == anchor:
            continue
        plus = sum(i in s and labels[j][s.index(i)] == 1
                   for j, s in enumerate(selectors))
        minus = sum(i in s and labels[j][s.index(i)] == -1
                    for j, s in enumerate(selectors))
        degree = plus + minus
        if degree:
            out += Fraction(2 * plus * minus, count * degree)
    return out


def brute_anchor_cost_a6(kappa: int = 20):
    """Brute-force the paired-orientation reduced game for A6,m=3,anchor=0."""
    a = A6
    n, m, anchor = 6, 3, 0
    q = qnorm(a)
    allowance = floor(((m / n) ** 1.5 - (m / n) ** 2) * q + 1e-12)
    selectors = [s for s in combinations(range(n), m) if anchor in s]
    choices = []
    for s in selectors:
        # Quotient duplicate projective labels after anchor orientation.
        unique = {}
        for row in favorable_rows(a, s, allowance):
            label = anchor_orient(row["label"], s, anchor)
            key = tuple(map(int, label))
            unique[key] = min(unique.get(key, q + 1), q - row["width"])
        choices.append([(np.asarray(key, dtype=np.int64), deficit)
                        for key, deficit in sorted(unique.items())])
    total = 1
    for c in choices:
        total *= len(c)
    best = None
    best_payload = None
    for assignment in product(*choices):
        labels = [x[0] for x in assignment]
        mean_deficit = Fraction(sum(x[1] for x in assignment), len(selectors))
        conflict = corrected_conflict(selectors, labels, anchor)
        cost = mean_deficit + kappa * conflict
        if best is None or cost < best:
            best = cost
            best_payload = (mean_deficit, conflict)
    assert best is not None
    return {
        "assignment_count": total,
        "best_cost": best,
        "mean_width_deficit": best_payload[0],
        "conflict": best_payload[1],
        "target_budget": Fraction(kappa, 2),
    }


def a6_all_subfamilies(kappa: int = 20):
    """Exact full game for every anchored-selector subfamily on A6,m=3.

    Each selector has one favorable projective label.  Its two response
    orientations have opposite edge words, and their scores can be balanced
    inside the fractional cube, so the one-block hard-favorable value is 8.
    Consequently W_G^F(kappa)=8-kappa*C_G exactly (row-free).
    """
    a = A6
    n, m, anchor = 6, 3, 0
    q = qnorm(a)
    allowance = 1
    selectors = [s for s in combinations(range(n), m) if anchor in s]
    labels = []
    for s in selectors:
        rows = favorable_rows(a, s, allowance)
        assert len(rows) == 1
        row = rows[0]
        label = anchor_orient(row["label"], s, anchor)
        assert row["width"] == 8
        # Verify directly that a cube point balances the two orientations.
        # Their intercept difference is at most twice the number of edges.
        assert abs(row["zplus"] - row["zminus"]) <= 4 * (m * (m - 1) // 2)
        labels.append(label)

    by_size = {}
    target_passing = []
    for mask in range(1, 1 << len(selectors)):
        idx = [j for j in range(len(selectors)) if mask >> j & 1]
        ss = [selectors[j] for j in idx]
        yy = [labels[j] for j in idx]
        conflict = corrected_conflict(ss, yy, anchor)
        value = Fraction(8) - kappa * conflict
        target = Fraction(q) - Fraction(kappa, 2)
        passes = value >= target
        size = len(idx)
        row = by_size.setdefault(size, {"min_conflict": conflict,
                                        "max_conflict": conflict,
                                        "pass_count": 0,
                                        "count": 0})
        row["min_conflict"] = min(row["min_conflict"], conflict)
        row["max_conflict"] = max(row["max_conflict"], conflict)
        row["pass_count"] += int(passes)
        row["count"] += 1
        if passes:
            target_passing.append(mask)
    return {
        "formula": "W=8-kappa*C",
        "target_conflict_threshold": Fraction(2, 5),
        "largest_passing_size": max(bin(mask).count("1") for mask in target_passing),
        "by_size": by_size,
    }


def anchored_exact_degree(a: np.ndarray, m: int, allowance: Optional[int] = None):
    """Maximum exact-favorable degree in each singleton-anchored slice."""
    n = len(a)
    q = qnorm(a)
    if allowance is None:
        allowance = floor(((m / n) ** 1.5 - (m / n) ** 2) * q + 1e-12)
    full = [x_of(mask, n) for mask in range(1 << (n - 1))]
    out = {}
    for anchor in range(n):
        selectors = [s for s in combinations(range(n), m) if anchor in s]
        favorable = {}
        for s in selectors:
            favorable[s] = {tuple(map(int, anchor_orient(row["label"], s, anchor)))
                             for row in favorable_rows(a, s, allowance)}
        best = 0
        for z0 in full:
            # Reorient the full projective word to +1 at the chosen anchor.
            z = z0 * int(z0[anchor])
            degree = sum(tuple(map(int, z[list(s)])) in favorable[s]
                         for s in selectors)
            best = max(best, degree)
        out[anchor] = Fraction(best, len(selectors))
    return out


def full_slice_game_upper_bounds(a: np.ndarray, m: int, allowance: Optional[int] = None,
                                 require_failure: bool = True):
    """Exact rational upper bounds for the full singleton-anchored games.

    For each block choose an optimizer of its unpenalized hard-favorable
    response game.  Every assignment then loses the mean one-block gap.  Its
    conflict is at least one minus the best exact center degree.
    """
    n = len(a)
    q = qnorm(a)
    kappa = 4 * (n - 1)
    if allowance is None:
        allowance = floor(((m / n) ** 1.5 - (m / n) ** 2) * q + 1e-12)
    degrees = anchored_exact_degree(a, m, allowance)
    out = {}
    for anchor in range(n):
        selectors = [s for s in combinations(range(n), m) if anchor in s]
        gaps = [q - hard_favorable_value(a, s, allowance) for s in selectors]
        mean_gap = Fraction(sum(gaps), len(gaps))
        conflict_lb = 1 - degrees[anchor]
        value_ub = Fraction(q) - mean_gap - kappa * conflict_lb
        target = Fraction(q) - Fraction(kappa, 2)
        out[anchor] = {
            "mean_gap": mean_gap,
            "conflict_lb": conflict_lb,
            "W_upper": value_ub,
            "target": target,
            "fails": value_ub < target,
        }
        if require_failure:
            assert value_ub < target
    return out


def slack_monotonicity_audits():
    """Show that enlarging favorable lists weakens both exact obstructions."""
    out = {}
    for a, name, m, allowances in (
        (A6, "A6", 3, (0, 4, 8)),
        (A8, "A8", 5, (0, 4, 8, 12)),
        (A9, "A9", 6, (0, 4, 8, 12, 16)),
    ):
        rows = []
        prior_gap, prior_degree = None, None
        for allowance in allowances:
            bounds = full_slice_game_upper_bounds(
                a, m, allowance=allowance, require_failure=False
            )
            max_gap = max(row["mean_gap"] for row in bounds.values())
            min_degree = min(1 - row["conflict_lb"] for row in bounds.values())
            if prior_gap is not None:
                assert max_gap <= prior_gap
                assert min_degree >= prior_degree
            prior_gap, prior_degree = max_gap, min_degree
            rows.append((allowance, max_gap, min_degree,
                         all(row["fails"] for row in bounds.values())))
        out[name] = rows
    return out


def coherent_width_certificate(a: np.ndarray, m: int, allowance: int):
    """Find an exact conflict-zero width certificate on a full anchored slice."""
    n = len(a)
    q = qnorm(a)
    kappa = 4 * (n - 1)
    full = [x_of(mask, n) for mask in range(1 << (n - 1))]
    out = {}
    for anchor in range(n):
        selectors = [s for s in combinations(range(n), m) if anchor in s]
        tables = {}
        for s in selectors:
            table = {}
            for row in favorable_rows(a, s, allowance):
                label = tuple(map(int, anchor_orient(row["label"], s, anchor)))
                table[label] = min(table.get(label, q + 1), q - row["width"])
            tables[s] = table
        best = None
        for z0 in full:
            z = z0 * int(z0[anchor])
            labels = [tuple(map(int, z[list(s)])) for s in selectors]
            if all(label in tables[s] for label, s in zip(labels, selectors)):
                cost = Fraction(sum(tables[s][label]
                                    for label, s in zip(labels, selectors)),
                                len(selectors))
                best = cost if best is None else min(best, cost)
        assert best is not None
        assert best <= Fraction(kappa, 2)
        out[anchor] = {"mean_width_deficit": best,
                       "target_budget": Fraction(kappa, 2)}
    return out


def main():
    for args in ((A6, "A6", 3), (A8, "A8", 5), (A9, "A9", 6)):
        print(audit_widths(*args))
        print(args[1], "anchored exact degrees", anchored_exact_degree(args[0], args[2]))
        print(args[1], "full-slice W upper bounds", full_slice_game_upper_bounds(args[0], args[2]))
    print("A6 anchor paired reduction", brute_anchor_cost_a6())
    print("A6 all anchored subfamilies", a6_all_subfamilies())
    print("slack monotonicity", slack_monotonicity_audits())
    print("slack coherent A6", coherent_width_certificate(A6, 3, 4))
    print("slack coherent A8", coherent_width_certificate(A8, 5, 8))
    print("slack coherent A9", coherent_width_certificate(A9, 6, 16))
    print("simultaneous_replacement_r36_check: PASS")


if __name__ == "__main__":
    main()
