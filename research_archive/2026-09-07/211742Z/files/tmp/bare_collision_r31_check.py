#!/usr/bin/env python3
"""Exact small checks for the Wave 31 bare-collision memo."""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product
from math import ceil, comb, e, log

import numpy as np
from scipy.optimize import linprog


def generic_moment_check() -> None:
    # A deliberately overlapping finite incidence system.
    sets = [
        frozenset((0, 1, 2)),
        frozenset((1, 2, 3)),
        frozenset((0, 3, 4)),
        frozenset((2, 4)),
    ]
    weights = [Fraction(1, 10), Fraction(2, 10), Fraction(3, 10),
               Fraction(1, 10), Fraction(3, 10)]
    M = len(sets)
    r, s = 4, 3
    p = [sum(weights[x] for x in H) for H in sets]
    S_s = sum(x**s for x in p)
    S_r = sum(x**r for x in p)

    event = Fraction(0)
    common = Fraction(0)
    ez2 = Fraction(0)
    for xs in product(range(len(weights)), repeat=r):
        mass = np.prod([weights[x] for x in xs])
        hits = [sum(x in H for x in xs) for H in sets]
        if max(hits) >= s:
            event += mass
        z = sum(all(x in H for x in xs) for H in sets)
        if z:
            common += mass
        ez2 += mass * z * z

    assert S_s / M <= event <= comb(r, s) * S_s
    assert common <= S_r
    assert common * ez2 >= S_r * S_r
    print("PASS generic exact sandwich and Paley--Zygmund", event, common)


def t_subset_wall_check(N: int = 7, t: int = 3, r: int = 6, s: int = 3) -> None:
    candidates = list(combinations(range(N), t))
    candidate_sets = [frozenset(B) for B in candidates]
    M = len(candidates)

    event_rs_count = 0
    event_rr_count = 0
    support_formula_count = sum(
        comb(N, u) * onto_count(r, u) for u in range(1, t + 1)
    )
    for xs in product(range(N), repeat=r):
        if len(set(xs)) <= t:
            event_rr_count += 1
        if any(sum(x in B for x in xs) >= s for B in candidate_sets):
            event_rs_count += 1
    assert event_rr_count == support_formula_count

    p = Fraction(t, N)
    S_s = M * p**s
    q_rs = Fraction(event_rs_count, N**r)
    q_rr = Fraction(event_rr_count, N**r)
    assert S_s / M <= q_rs <= comb(r, s) * S_s
    assert q_rr <= M * p**r

    # Numerical LP audit of delta*=t/N.
    H = np.zeros((M, N))
    for j, B in enumerate(candidate_sets):
        H[j, list(B)] = 1.0
    # Variables w_1,...,w_N,v; minimize v with H w <= v.
    Aub = np.hstack((H, -np.ones((M, 1))))
    Aeq = np.zeros((1, N + 1))
    Aeq[0, :N] = 1.0
    objective = np.zeros(N + 1)
    objective[-1] = 1.0
    sol = linprog(
        objective,
        A_ub=Aub,
        b_ub=np.zeros(M),
        A_eq=Aeq,
        b_eq=np.ones(1),
        bounds=[(0.0, None)] * N + [(0.0, 1.0)],
        method="highs",
    )
    assert sol.success
    assert abs(sol.fun - t / N) < 1e-10
    assert Fraction(comb(N - 1, t - 1), M) == p
    print("PASS t-subset wall and exact game value", {
        "N": N,
        "t": t,
        "r": r,
        "s": s,
        "J_rr": str(q_rr),
        "J_rs": str(q_rs),
        "delta": str(p),
    })


def growing_group_algebra_check() -> None:
    # Audit the constants in (R31.2b) with log M at its allowed maximum.
    r, T, L0, K = 53, 5, 7.0, 2.0
    s = ceil(r / T)
    log_M = r * L0
    log_tau_bound_from_sandwich = (
        log_M + log(comb(r, s)) + K * r * L0
    ) / s
    claimed = (K + 1) * T * L0 + log(e * T)
    assert r / s <= T
    assert log_tau_bound_from_sandwich <= claimed + 1e-12
    print("PASS growing-group constants", {
        "r": r, "s": s, "T": T,
        "direct": log_tau_bound_from_sandwich,
        "claimed": claimed,
    })


def onto_count(r: int, u: int) -> int:
    # u! S(r,u), by inclusion--exclusion.
    return sum((-1) ** j * comb(u, j) * (u - j) ** r for j in range(u + 1))


if __name__ == "__main__":
    generic_moment_check()
    t_subset_wall_check()
    growing_group_algebra_check()
    print("PASS all Wave 31 bare-collision checks")
