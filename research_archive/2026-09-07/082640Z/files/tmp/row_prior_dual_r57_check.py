#!/usr/bin/env python3
"""Finite checks for the Wave 57 row-penalized favorable-incidence dual.

All signing energies and row costs are enumerated exactly.  The threshold
contains (m/n)^(3/2), so favorable decisions use floating point; we report and
assert a positive margin from every decision boundary.
"""

from __future__ import annotations

import itertools
import math
import sys

import numpy as np
from scipy.optimize import linprog

sys.path.insert(0, "/home/math/quadra/tmp")
from envelope_block_cover_r27 import A6, A8, A9, projective_spins


def qnorm(a: np.ndarray) -> int:
    return max(abs(int(x @ a @ x)) for x in projective_spins(len(a)))


def columns(a: np.ndarray, m: int, t: float = 0.0):
    """Return (coverage u, surprise h, row, deficit) for every oriented cut."""
    n = len(a)
    q = qnorm(a)
    selectors = list(itertools.combinations(range(n), m))
    q_s = []
    for s in selectors:
        b = a[np.ix_(s, s)]
        q_s.append(qnorm(b))
    p2 = m * (m - 1) / (n * (n - 1))
    b_allow = ((m / n) ** 1.5 - p2) * q
    ans = []
    margin = math.inf
    for x in projective_spins(n):
        raw = int(x @ a @ x)
        row = int((a @ x) @ (a @ x))
        for sigma in (-1, 1):
            deficit = q - sigma * raw
            count = 0
            for s, qs in zip(selectors, q_s):
                local = sigma * int(x[list(s)] @ a[np.ix_(s, s)] @ x[list(s)])
                slack = b_allow + t + p2 * deficit - (qs - local)
                margin = min(margin, abs(slack))
                count += slack >= 0
            u = count / len(selectors)
            if count:
                ans.append((u, -math.log(u), float(row), float(deficit)))
    assert margin > 1e-9
    return np.asarray(ans), margin


def pythagorean_check() -> None:
    """Check K=E h+E posterior KL on a nonuniform abstract incidence."""
    u0 = np.asarray([0.2, 0.5, 0.8])
    # Uniform selector has ten atoms; fibres are initial subsets of sizes 2,5,8.
    fibres = [np.arange(k) for k in (2, 5, 8)]
    pi = np.asarray([0.2, 0.3, 0.5])
    conditionals = []
    for j, f in enumerate(fibres):
        raw = np.arange(1, len(f) + 1, dtype=float) ** (j + 1)
        p = np.zeros(10)
        p[f] = raw / raw.sum()
        conditionals.append(p)
    joint = np.stack([pi[j] * conditionals[j] for j in range(3)], axis=1)
    ps = joint.sum(axis=1)
    pd = joint.sum(axis=0)
    mutual = sum(
        joint[s, d] * math.log(joint[s, d] / (ps[s] * pd[d]))
        for s in range(10) for d in range(3) if joint[s, d] > 0
    )
    selector_kl = sum(x * math.log(x / 0.1) for x in ps if x > 0)
    rhs = 0.0
    for d, f in enumerate(fibres):
        posterior_kl = sum(
            conditionals[d][s]
            * math.log(conditionals[d][s] / (1.0 / len(f)))
            for s in f if conditionals[d][s] > 0
        )
        rhs += pi[d] * (-math.log(u0[d]) + posterior_kl)
    assert abs(mutual + selector_kl - rhs) < 1e-12


def frontier_dual_check(data: np.ndarray, h_budget: float) -> dict:
    """Primal/dual for min E row subject to E surprise <= h_budget."""
    h = data[:, 1]
    row = data[:, 2]
    k = len(data)
    primal = linprog(
        row,
        A_ub=h[None, :],
        b_ub=np.asarray([h_budget]),
        A_eq=np.ones((1, k)),
        b_eq=np.ones(1),
        bounds=(0, None),
        method="highs",
    )
    assert primal.success
    # max_z z-lambda*H subject to z<=row_d+lambda*h_d, lambda>=0.
    dual = linprog(
        np.asarray([h_budget, -1.0]),
        A_ub=np.column_stack((-h, np.ones(k))),
        b_ub=row,
        bounds=((0, None), (None, None)),
        method="highs",
    )
    assert dual.success
    dual_value = -dual.fun
    assert abs(primal.fun - dual_value) < 1e-8
    support = np.flatnonzero(primal.x > 1e-8)
    assert len(support) <= 2
    return {
        "h_budget": h_budget,
        "min_mean_row": primal.fun,
        "support": [
            (float(primal.x[i]), float(h[i]), float(row[i])) for i in support
        ],
        "dual_lambda": float(dual.x[0]),
    }


def selected_prior_lp_check(data: np.ndarray, row_budget: float) -> dict:
    """Max selected-prior incidence with a captured mean-row constraint."""
    u = data[:, 0]
    row = data[:, 2]
    k = len(data)
    # q_d=nu_d*u_d.  Then sum q/u=1 and captured row is <= row_budget.
    primal = linprog(
        -np.ones(k),
        A_ub=(row - row_budget)[None, :],
        b_ub=np.zeros(1),
        A_eq=(1.0 / u)[None, :],
        b_eq=np.ones(1),
        bounds=(0, None),
        method="highs",
    )
    assert primal.success
    z = -primal.fun
    q = primal.x
    nu = q / u
    assert abs(nu.sum() - 1.0) < 1e-9
    captured = q / z
    assert captured @ row <= row_budget + 1e-8
    # Direct dual: min_y y, y/u_d+lambda(row_d-row_budget)>=1.
    dual = linprog(
        np.asarray([1.0, 0.0]),
        A_ub=np.column_stack((-1.0 / u, -(row - row_budget))),
        b_ub=-np.ones(k),
        bounds=((None, None), (0, None)),
        method="highs",
    )
    assert dual.success
    assert abs(z - dual.fun) < 1e-8
    support = np.flatnonzero(q > 1e-9)
    assert len(support) <= 2
    # Posterior parametrization: Z^{-1}=E_captured(1/u).
    assert abs(1.0 / z - captured @ (1.0 / u)) < 1e-9
    return {
        "row_budget": row_budget,
        "max_incidence": z,
        "support": [
            (float(captured[i]), float(-math.log(u[i])), float(row[i]))
            for i in support
        ],
        "dual_row_price": float(dual.x[1]),
    }


def fixed_prior_kl_check(data: np.ndarray) -> None:
    """Check the Gibbs/KL identity for a fixed nonuniform selected prior."""
    k = len(data)
    raw = np.arange(1, k + 1, dtype=float)
    nu = raw / raw.sum()
    u, row = data[:, 0], data[:, 2]
    lam = 0.017
    z = float(np.sum(nu * u * np.exp(-lam * row)))
    # Optimizing joint law: uniform inside each fibre and tilted output.
    pi = nu * u * np.exp(-lam * row) / z
    kcost = float(pi @ (-np.log(u)))
    outkl = float(sum(x * math.log(x / y) for x, y in zip(pi, nu) if x))
    objective = kcost + outkl + lam * float(pi @ row)
    assert abs(objective + math.log(z)) < 2e-12


def sector_bound_check(a: np.ndarray) -> None:
    """Check R<=2(n-1) sum_i (signed field_i)_+ when energy is nonnegative."""
    n = len(a)
    for x in projective_spins(n):
        raw = int(x @ a @ x)
        for sigma in (-1, 1):
            fields = sigma * x * (a @ x)
            energy = sigma * raw
            assert int(fields.sum()) == energy
            row = int(fields @ fields)
            if energy >= 0:
                positive_l1 = int(np.maximum(fields, 0).sum())
                assert row <= 2 * (n - 1) * positive_l1


def a8_orientation_wall() -> None:
    """Exact wall: the largest A8 bare fibres are opposite-oriented grounds."""
    a, m = A8, 5
    n = len(a)
    q = qnorm(a)
    selectors = list(itertools.combinations(range(n), m))
    q_s = [qnorm(a[np.ix_(s, s)]) for s in selectors]
    p2 = m * (m - 1) / (n * (n - 1))
    b_allow = ((m / n) ** 1.5 - p2) * q
    records = []
    for x in projective_spins(n):
        raw = int(x @ a @ x)
        row = int((a @ x) @ (a @ x))
        counts = {}
        for sigma in (-1, 1):
            deficit = q - sigma * raw
            count = sum(
                qs - sigma * int(x[list(s)] @ a[np.ix_(s, s)] @ x[list(s)])
                <= b_allow + p2 * deficit + 1e-12
                for s, qs in zip(selectors, q_s)
            )
            fields = sigma * x * (a @ x)
            counts[sigma] = count
            records.append((count, row, deficit, fields, x.copy(), sigma))
    best = max(z[0] for z in records)
    maximizers = [z for z in records if z[0] == best]
    assert (q, best, len(selectors), len(maximizers)) == (20, 40, 56, 8)
    assert all(z[1] == 64 and z[2] == 2 * q for z in maximizers)
    assert all(int(np.maximum(z[3], 0).sum()) == 0 for z in maximizers)
    # Reversing only the orientation makes each word an exact parent ground,
    # keeps its row, and lowers rather than raises its bare coverage.
    for _, row, _, _, x, sigma in maximizers:
        mate = next(
            z for z in records
            if z[5] == -sigma and np.array_equal(z[4], x)
        )
        assert mate[1] == row and mate[2] == 0 and mate[0] == 30


def hub_avoidance_scale_check() -> None:
    """Numerically audit the sharp scalable entropy/row exponent obstruction."""
    # c=1/8, k=n^(1/4-c)*log n.  Then k=o(H), while k*n^2 is
    # omega(n^(9/4-c)).  The exact log mass is evaluated by a stable sum.
    c = 1.0 / 8.0
    ratios = []
    for n in (10_000, 100_000, 1_000_000):
        p = 0.6
        m = int(p * n)
        k = max(1, int(n ** (0.25 - c) * math.log(n)))
        log_inverse_mass = sum(
            math.log((n - j) / (n - m - j)) for j in range(k)
        )
        h_scale = n ** (0.75 - c)
        row_ratio = (k * n * n) / (n ** (2.25 - c))
        ratios.append((log_inverse_mass / h_scale, row_ratio))
    assert ratios[-1][0] < ratios[0][0]
    assert ratios[-1][1] > ratios[0][1]


def audit(name: str, a: np.ndarray, m: int) -> None:
    data, margin = columns(a, m)
    fixed_prior_kl_check(data)
    sector_bound_check(a)
    hmin, hmax = float(data[:, 1].min()), float(data[:, 1].max())
    rmin, rmax = float(data[:, 2].min()), float(data[:, 2].max())
    budgets = sorted(set((hmin, 0.5 * (hmin + hmax), hmax)))
    frontiers = [frontier_dual_check(data, b) for b in budgets]
    row_budgets = sorted(set((rmin, 0.5 * (rmin + rmax), rmax)))
    selected = [selected_prior_lp_check(data, b) for b in row_budgets]
    print({
        "name": name,
        "m": m,
        "decision_margin": margin,
        "column_count": len(data),
        "h_range": (hmin, hmax),
        "row_range": (rmin, rmax),
        "frontiers": frontiers,
        "selected_prior": selected,
    })


def main() -> None:
    pythagorean_check()
    a8_orientation_wall()
    hub_avoidance_scale_check()
    for a, name, m in ((A6, "A6", 4), (A8, "A8", 5), (A9, "A9", 6)):
        audit(name, a, m)
    print("PASS row-incidence projection, LP duality, and fixed-prior KL duality")


if __name__ == "__main__":
    main()
