#!/usr/bin/env python3
"""Exact finite audits for Wave 37 completion-width replacement candidates."""

from fractions import Fraction
from itertools import combinations

import numpy as np

from coherent_signatures_r29_search import A6, A8, A9
from simultaneous_replacement_r36_check import qnorm, response_table


QMIN = {3: 6, 5: 8, 6: 10}


def selector_data(a: np.ndarray, m: int):
    q = qnorm(a)
    qm = QMIN[m]
    rows = []
    for s in combinations(range(len(a)), m):
        table = response_table(a, s)
        qs = max(abs(r["energy"]) for r in table)
        grounds = [r for r in table if abs(r["energy"]) == qs]
        best_b = min(q - r["width"] for r in grounds)
        all_b = [q - r["width"] for r in grounds]
        outside = tuple(i for i in range(len(a)) if i not in s)
        qo = 0 if not outside else qnorm(a[np.ix_(outside, outside)])
        rows.append((s, qs, qo, best_b, min(all_b), max(all_b)))
    return q, qm, rows


def audit(a: np.ndarray, name: str, m: int):
    q, qm, rows = selector_data(a, m)
    # Candidate pointwise and averaged inequalities.
    point_gap = [(s, b - (q - qm)) for s, qs, qo, b, _, _ in rows]
    excesses = [g for _, g in point_gap]
    # Exact zero-block/interval-center bound proved in the memo:
    # min over all labels b <= qm + Q(A[S^c]).  Check it by enumeration.
    zero_bounds = []
    for s in combinations(range(len(a)), m):
        table = response_table(a, s)
        best_all = min(q - r["width"] for r in table)
        outside = tuple(i for i in range(len(a)) if i not in s)
        qo = 0 if not outside else qnorm(a[np.ix_(outside, outside)])
        assert best_all <= qm + qo
        zero_bounds.append((best_all, qm + qo))

    # Test stronger favorable variants of the zero-block inequality.
    fav_zero_fail = [(s, b, qm + qo) for s, qs, qo, b, _, _ in rows
                     if b > qm + qo]
    excess_adjusted = [b - (q - qm + (qs - qm))
                       for s, qs, qo, b, _, _ in rows]

    # Two-sign self-exposure dichotomy, checked for every exact child ground.
    exposure = []
    for s in combinations(range(len(a)), m):
        table = response_table(a, s)
        qs = max(abs(r["energy"]) for r in table)
        for row in table:
            if abs(row["energy"]) != qs:
                continue
            b = q - row["width"]
            center = Fraction(row["zplus"] - row["zminus"], 2)
            eps = max(Fraction(0), b - abs(center + qm), b - abs(center - qm))
            delta = qs - qm
            assert eps >= delta or 2 * b <= q - qm + eps
            exposure.append((eps, delta, b))
    return {
        "name": name,
        "q": q,
        "qm": qm,
        "best_ground_b_hist": {v: sum(r[3] == v for r in rows)
                               for v in sorted({r[3] for r in rows})},
        "max_point_excess_over_q_minus_qm": max(excesses),
        "mean_point_excess": Fraction(sum(excesses), len(excesses)),
        "point_fail_count": sum(g > 0 for g in excesses),
        "adjusted_qs_fail_count": sum(g > 0 for g in excess_adjusted),
        "adjusted_qs_max_excess": max(excess_adjusted),
        "self_exposure_epsilon_range": (min(x[0] for x in exposure),
                                         max(x[0] for x in exposure)),
        "self_exposure_migration_count": sum(eps >= delta and delta > 0
                                              for eps, delta, b in exposure),
        "favorable_zero_bound_fail_count": len(fav_zero_fail),
        "first_favorable_zero_failure": fav_zero_fail[:1],
        "zero_bound_max_slack": max(rhs - lhs for lhs, rhs in zero_bounds),
    }


def main():
    for args in ((A6, "A6", 3), (A8, "A8", 5), (A9, "A9", 6)):
        print(audit(*args))
    print("width_replacement_r37_check: PASS")


if __name__ == "__main__":
    main()
