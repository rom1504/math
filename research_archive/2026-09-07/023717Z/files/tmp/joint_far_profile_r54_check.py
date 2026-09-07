#!/usr/bin/env python3
"""Finite audit of the Wave 54 orientation-pair/profile reduction.

This is diagnostic only: H=min(3,|T|), t=0, and b_H=0 are finite
normalizations, not the asymptotic far-tail parameters.
"""

from __future__ import annotations

import itertools
import math
import sys

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from cross_block_dual_r48_check import sample_order_ten
from envelope_block_cover_r27 import A6, A8, A9, projective_spins
from far_tail_kfunctional_r52_check import exact_kfunctional, qnorm


def huber_kfunctional(x: np.ndarray, hfun: int) -> float:
    """Evaluate the exact one-dimensional Huber variational formula."""

    b = np.abs(np.asarray(x, dtype=float))
    support = int(np.count_nonzero(b))
    if support <= hfun:
        return float(np.sum(b))
    lo = 0.0
    hi = float(np.max(b))
    while float(np.sum(np.minimum(1.0, b / hi) ** 2)) > hfun:
        hi *= 2.0
    for _ in range(120):
        theta = (lo + hi) / 2.0
        if float(np.sum(np.minimum(1.0, b / theta) ** 2)) > hfun:
            lo = theta
        else:
            hi = theta
    theta = hi
    phi = np.where(
        b <= theta,
        b * b / (2 * theta),
        b - theta / 2,
    )
    return hfun * theta / 2 + float(np.sum(phi))


def audit_huber_formula() -> int:
    rng = np.random.default_rng(5401)
    checked = 0
    for k in range(1, 18):
        for hfun in range(1, k + 2):
            for _ in range(20):
                x = rng.integers(-12, 13, size=k).astype(float)
                exact = exact_kfunctional(x, hfun)
                huber = huber_kfunctional(x, hfun)
                assert abs(exact - huber) < 1e-8
                checked += 1
    return checked


def audit_instance(a: np.ndarray, name: str, m: int) -> dict[str, int | float | str]:
    n = len(a)
    q = qnorm(a)
    p = m / n
    p2 = m * (m - 1) / (n * (n - 1))
    hfun = min(3, n - m)

    oriented = 0
    negative = 0
    partner_nonnegative = 0
    both_negative = 0
    both_negative_profile_success = 0
    both_negative_profile_failure = 0
    pair_identity_checks = 0
    pair_window_checks = 0
    min_ratio = math.inf

    for s in itertools.combinations(range(n), m):
        tset = tuple(i for i in range(n) if i not in s)
        child = a[np.ix_(s, s)]
        cross = a[np.ix_(tset, s)]
        h = qnorm(child) - p ** 1.5 * q
        for y0 in projective_spins(m):
            y = np.asarray(y0, dtype=np.int64)
            a0 = int(y @ child @ y)
            u = cross @ y
            kval = exact_kfunctional(u.astype(float), hfun)
            gs = {
                sigma: (1 - p2) * sigma * a0 - h
                for sigma in (-1, 1)
            }
            d = (1 - p2) * abs(a0)
            for far_cut in (0.5, 2.0, 5.0, 11.0):
                lower = p2 * far_cut
                upper = p2 * kval
                predicted = int(lower <= h - d <= upper) + int(
                    lower <= h + d <= upper
                )
                actual = 0
                for sigma_test in (-1, 1):
                    g_test = gs[sigma_test]
                    if g_test < 0:
                        r_test = -g_test / p2
                        actual += int(r_test >= far_cut and kval >= r_test)
                assert predicted == actual
                pair_window_checks += 1
            for sigma in (-1, 1):
                oriented += 1
                g = gs[sigma]
                if g >= 0:
                    continue
                negative += 1
                r = -g / p2
                gopp = gs[-sigma]
                assert abs(gopp - (-g - 2 * h)) < 1e-10
                assert abs(gopp - (p2 * r - 2 * h)) < 1e-10
                pair_identity_checks += 1
                if gopp >= 0:
                    partner_nonnegative += 1
                else:
                    both_negative += 1
                    if kval + 1e-10 >= r:
                        both_negative_profile_success += 1
                    else:
                        both_negative_profile_failure += 1
                        if kval > 0:
                            min_ratio = min(min_ratio, kval / r)
                        else:
                            min_ratio = 0.0

    assert negative == partner_nonnegative + both_negative
    return {
        "name": name,
        "n": n,
        "m": m,
        "q": q,
        "H": hfun,
        "oriented": oriented,
        "negative": negative,
        "partner_nonnegative": partner_nonnegative,
        "both_negative": both_negative,
        "both_negative_profile_success": both_negative_profile_success,
        "both_negative_profile_failure": both_negative_profile_failure,
        "min_failure_K_over_r": min_ratio,
        "pair_identity_checks": pair_identity_checks,
        "pair_window_checks": pair_window_checks,
    }


def main() -> None:
    huber_checks = audit_huber_formula()
    print(f"Huber/exact K-functional checks={huber_checks}")
    rows: list[dict[str, int | float | str]] = []
    cases = []
    for a, name in ((A6, "A6"), (A8, "A8"), (A9, "A9"),
                    (sample_order_ten(0), "A10_seed0")):
        n = len(a)
        for m in range(n // 2 + 1, n):
            cases.append((a, name, m))
    for a, name, m in cases:
        row = audit_instance(a, name, m)
        rows.append(row)
        print(
            f"{name} m={m} q={row['q']} H={row['H']} "
            f"negative={row['negative']}/{row['oriented']} "
            f"partner>=0={row['partner_nonnegative']} "
            f"both<0={row['both_negative']} "
            f"both-branch K>=r={row['both_negative_profile_success']} "
            f"K<r={row['both_negative_profile_failure']} "
            f"min(K/r|fail)={row['min_failure_K_over_r']:.6g}"
        )
    assert sum(int(r["both_negative"]) for r in rows) > 0
    assert sum(int(r["both_negative_profile_failure"]) for r in rows) > 0
    print(
        "TOTAL negative/partner>=0/both<0/K-success/K-failure=",
        sum(int(r["negative"]) for r in rows),
        sum(int(r["partner_nonnegative"]) for r in rows),
        sum(int(r["both_negative"]) for r in rows),
        sum(int(r["both_negative_profile_success"]) for r in rows),
        sum(int(r["both_negative_profile_failure"]) for r in rows),
    )
    print(
        "TOTAL pair identity checks=",
        sum(int(r["pair_identity_checks"]) for r in rows),
    )
    print(
        "TOTAL pair-window checks=",
        sum(int(r["pair_window_checks"]) for r in rows),
    )
    print("PASS Wave 54 orientation-pair finite audit")


if __name__ == "__main__":
    main()
