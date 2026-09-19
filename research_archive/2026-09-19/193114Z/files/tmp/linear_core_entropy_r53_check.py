#!/usr/bin/env python3
"""Deterministic checks for the Wave 53 linear-core collision audit."""

from __future__ import annotations

import itertools
import math
import sys
from collections import Counter
from fractions import Fraction

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from box_triple_r50_check import A as A10
from coarea_alignment_r47_check import exact_boundary_data


def falling(x: int, j: int) -> int:
    ans = 1
    for i in range(j):
        ans *= x - i
    return ans


def eigenvalue(n: int, m: int, ell: int, j: int) -> Fraction:
    if ell < j or n - m < j:
        return Fraction(0)
    return Fraction(
        falling(ell, j) * falling(n - m, j),
        falling(m, j) * falling(n - ell, j),
    )


def core_statistics(n: int, m: int, family, ell: int):
    states = [frozenset(s) for s in family]
    r = len(states)
    assert r
    loads = Counter(
        core
        for s in states
        for core in itertools.combinations(sorted(s), ell)
    )
    intersections = Counter(len(s & t) for s in states for t in states)
    d = math.comb(m, ell) * math.comb(n - ell, m - ell)
    b = (
        math.comb(m - 2, ell - 2)
        * math.comb(n - ell - 2, m - ell)
    )
    diagonal = math.comb(m, ell)
    off_collision = sum(c * (c - 1) for c in loads.values())
    off_histogram = sum(
        count * math.comb(j, ell)
        for j, count in intersections.items()
        if j < m
    )
    overlap_partner_count = sum(
        count for j, count in intersections.items() if ell <= j < m
    )
    assert off_collision == off_histogram
    delta = Fraction(off_collision, r)
    partners = Fraction(overlap_partner_count, r)
    assert delta <= math.comb(m, ell) * partners
    p_retention = Fraction(sum(c * c for c in loads.values()), r * d)
    h = Fraction(diagonal, d)
    lambda2 = eigenvalue(n, m, ell, 2)
    assert h == Fraction(1, math.comb(n - ell, m - ell))
    assert lambda2 == Fraction(b, d)
    assert p_retention == h + delta / d
    assert p_retention - lambda2 == Fraction(
        off_collision - r * (b - diagonal), r * d
    )

    q_denom = r * math.comb(m, ell)
    collision_q = sum(Fraction(c, q_denom) ** 2 for c in loads.values())
    k_eff = Fraction(1, 1) / collision_q
    a = Fraction(r, math.comb(n, m))
    assert p_retention == a * Fraction(math.comb(n, ell), 1) / k_eff
    return {
        "r": r,
        "a": a,
        "p": p_retention,
        "h": h,
        "lambda2": lambda2,
        "delta": delta,
        "partners": partners,
        "threshold": b - diagonal,
        "d": d,
        "k_eff": k_eff,
        "load_histogram": Counter(loads.values()),
        "intersection_histogram": intersections,
    }


def exhaustive_identity_check() -> None:
    n, m = 5, 3
    universe = [frozenset(s) for s in itertools.combinations(range(n), m)]
    checked = 0
    for mask in range(1, 1 << len(universe)):
        family = [universe[i] for i in range(len(universe)) if mask >> i & 1]
        for ell in (2,):
            core_statistics(n, m, family, ell)
            checked += 1
    print(f"exhaustive weighted-neighbor/entropy identities: {checked} families")


def a10_row_ten_check() -> None:
    n, m = 10, 6
    _, selectors, _, _, rows, _, _, favorable, degree = exact_boundary_data(A10, m)
    active = [i for i, az in enumerate(degree) if az > 0]
    cap = min(int(rows[i]) for i in active)
    low = [i for i in active if int(rows[i]) == cap]
    assert cap == 10 and len(low) == 2
    families = [tuple(np.flatnonzero(favorable[i])) for i in low]
    assert families[0] == families[1] and len(families[0]) == 5
    family = [frozenset(selectors[i]) for i in families[0]]
    print("A10,m=6,row=10 exact-ground family")
    for ell in range(2, m):
        data = core_statistics(n, m, family, ell)
        # Only ell=2 has positive excess.  Every genuinely larger audited core
        # is below the second eigenvalue despite the exact row-optimal witness.
        if ell == 2:
            assert data["p"] > data["lambda2"]
        else:
            assert data["p"] < data["lambda2"]
        print(
            f"  ell={ell}: a={data['a']} h={data['h']} "
            f"p={data['p']} lambda2={data['lambda2']} "
            f"p-lambda2={data['p']-data['lambda2']} "
            f"off_weighted_degree={data['delta']} "
            f"overlap_partners={data['partners']} "
            f"threshold={data['threshold']} K_eff={data['k_eff']}"
        )
    print(
        "  ordered intersections="
        + str(dict(sorted(core_statistics(n, m, family, 2)["intersection_histogram"].items())))
    )


def linear_scale_check() -> None:
    # Rational fixed-density sequence: p=3/5, alpha=3/10, beta=9/20.
    p = 3 / 5
    alpha = 3 / 10
    beta = 9 / 20
    lambda2_limit = (alpha * (1 - p) / (p * (1 - alpha))) ** 2
    entropy = lambda x: -x * math.log(x) - (1 - x) * math.log(1 - x)
    i_rho = entropy(alpha) - p * entropy(alpha / p)
    i_gamma = p * entropy(alpha / p) - beta * entropy(alpha / beta)
    assert i_rho > i_gamma > 0
    print(
        "linear limits p=3/5,alpha=3/10,beta=9/20: "
        f"lambda2->{lambda2_limit:.12f}, I_rho={i_rho:.12f}, "
        f"I_Gamma={i_gamma:.12f}, I_rho-I_Gamma={i_rho-i_gamma:.12f}"
    )
    for n in (100, 200, 400, 800):
        m, ell, s = 3 * n // 5, 3 * n // 10, 9 * n // 20
        rho = Fraction(math.comb(m, ell), math.comb(n, ell))
        gamma = Fraction(math.comb(s, ell), math.comb(m, ell))
        lam2 = eigenvalue(n, m, ell, 2)
        h = Fraction(1, math.comb(n - ell, m - ell))
        observed_i_rho = -math.log(float(rho)) / n
        observed_i_ratio = -math.log(float(rho / gamma)) / n
        assert h < lam2
        print(
            f"  n={n}: lambda2={float(lam2):.10f} "
            f"-log(rho)/n={observed_i_rho:.10f} "
            f"-log(rho/Gamma)/n={observed_i_ratio:.10f} "
            f"-log(h)/n={-math.log(float(h))/n:.10f}"
        )


def main() -> None:
    exhaustive_identity_check()
    a10_row_ten_check()
    linear_scale_check()
    print("PASS: linear-core weighted collision and finite obstruction checks")


if __name__ == "__main__":
    main()
