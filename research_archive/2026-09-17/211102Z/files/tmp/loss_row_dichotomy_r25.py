#!/usr/bin/env python3
"""Wave 25 audits for the loss-Laplace/coarea row-cost dichotomy."""

from __future__ import annotations

import itertools
import math
from collections import Counter

import numpy as np
from scipy.optimize import linprog

from check_response_dual_r16 import A8, A9, spins
from verify_compatible_replacement_r12 import A6


def ground_records(A: np.ndarray, m: int):
    n = len(A)
    X = spins(n)
    raw = np.einsum("bi,ij,bj->b", X, A, X)
    q = int(np.max(np.abs(raw)))
    grounds = []
    for sigma in (-1, 1):
        for x, value in zip(X, raw):
            if sigma * value != q:
                continue
            row = sigma * x * (x @ A)
            assert np.all(row >= 0)
            assert int(row.sum()) == q
            grounds.append((sigma, x.copy(), row.astype(int), int(row @ row)))

    selectors = list(itertools.combinations(range(n), m))
    child_q = {}
    for S in selectors:
        C = A[np.ix_(S, S)]
        Y = spins(m)
        values = np.einsum("bi,ij,bj->b", Y, C, Y)
        child_q[S] = int(np.max(np.abs(values)))

    records = []
    k = n - m
    for sigma, x, row, c in grounds:
        losses = []
        deletion_rows = []
        for S in selectors:
            T = tuple(i for i in range(n) if i not in S)
            score = sigma * int(x[list(S)] @ A[np.ix_(S, S)] @ x[list(S)])
            loss = child_q[S] - score
            assert loss >= 0
            R = int(row[list(T)].sum())
            b = int(
                sigma
                * sum(A[i, j] * x[i] * x[j] for i in T for j in S)
            )
            d = q - child_q[S]
            assert loss == R + b - d
            h = int(sigma * x[list(T)] @ A[np.ix_(T, T)] @ x[list(T)])
            assert 2 * R == d + h + loss
            losses.append(loss)
            deletion_rows.append(R)

        # Exact fixed-size sampling moments of the nonnegative ground rows.
        deletion_rows = np.asarray(deletion_rows, dtype=float)
        expected_R = k * q / n
        expected_R2 = (
            k * c / n
            + k * (k - 1) * (q * q - c) / (n * (n - 1))
        )
        assert math.isclose(float(deletion_rows.mean()), expected_R, abs_tol=1e-10)
        assert math.isclose(float((deletion_rows**2).mean()), expected_R2, abs_tol=1e-10)
        records.append({"sigma": sigma, "x": x, "rows": row, "cost": c, "losses": np.asarray(losses)})
    return q, selectors, records


def captured(records, threshold: float, prior: np.ndarray):
    u = np.array([(record["losses"] <= threshold).mean() for record in records])
    c = np.array([record["cost"] for record in records], dtype=float)
    Z = float(prior @ u)
    Ccap = float(prior @ (u * c) / Z) if Z > 0 else math.inf
    return Z, Ccap, u


def soft_lp(records, lam: float, budget: float):
    v = np.array([np.exp(-lam * record["losses"]).mean() for record in records])
    c = np.array([record["cost"] for record in records], dtype=float)
    # Include the dummy action (v,c)=(0,0).
    vv = np.r_[v, 0.0]
    cc = np.r_[c, 0.0]
    result = linprog(
        -vv,
        A_ub=np.array([vv * (cc - budget)]),
        b_ub=np.array([0.0]),
        A_eq=np.ones((1, len(vv))),
        b_eq=np.ones(1),
        bounds=(0, None),
        method="highs",
    )
    assert result.success
    prior = result.x[:-1]
    K = float(prior @ v)
    C_lap = float(prior @ (v * c) / K) if K > 0 else math.inf
    return K, C_lap, prior, v, c


def coarea_audit(records, lam: float, budget: float):
    K, C_lap, prior, v, costs = soft_lp(records, lam, budget)
    if K == 0:
        return {"K": 0.0, "support": 0}
    assert C_lap <= budget + 1e-8
    # General layer-cake constants: with a>2, discard at most K/a in the
    # upper t-tail and at most K/a on levels whose captured cost exceeds
    # a*C.  The remaining weighted level mass is K(1-2/a).
    factor = 4.0
    retained = 1 - 2 / factor
    threshold_cap = math.log(factor / K) / lam
    possible = sorted(
        {
            float(loss)
            for record in records
            for loss in record["losses"]
            if loss <= threshold_cap + 1e-10
        }
    )
    witnesses = []
    for threshold in possible:
        Z, Ccap, _ = captured(records, threshold, prior)
        if Z + 1e-10 >= retained * K and Ccap <= factor * budget + 1e-8:
            witnesses.append((threshold, Z, Ccap))
    assert witnesses
    chosen = min(witnesses, key=lambda item: item[0])

    # R_2 is integer valued.  Rounding the witnessed hard budget upward and
    # using the captured law gives an individual low-cost column.  If p is
    # its captured mass on {c <= C_int}, then
    #   C_int >= E c >= (1-p)(C_int+1),
    # so p >= 1/(C_int+1); the harmonic identity then forces one such
    # column to have u >= Z/(C_int+1).
    threshold, Z, Ccap = chosen
    hard_budget = int(math.ceil(factor * budget))
    _, _, u = captured(records, threshold, prior)
    low = costs <= hard_budget
    captured_weights = prior * u / Z
    low_mass = float(captured_weights[low].sum())
    assert Ccap <= hard_budget + 1e-8
    assert low_mass + 1e-10 >= 1 / (hard_budget + 1)
    assert float(np.max(u[low])) + 1e-10 >= Z / (hard_budget + 1)

    # Numerically verify the Laplace coarea identity ground by ground using
    # the discrete Stieltjes form E exp(-lambda L).
    for record, value in zip(records, v):
        direct = float(np.exp(-lam * record["losses"]).mean())
        assert math.isclose(value, direct, rel_tol=1e-12)

    return {
        "K": K,
        "laplace_captured_cost": C_lap,
        "support": int(np.sum(prior > 1e-9)),
        "threshold_cap": threshold_cap,
        "coarea_witness": chosen,
        "integer_extraction": {
            "budget": hard_budget,
            "low_captured_mass": low_mass,
            "coverage_lower_bound": Z / (hard_budget + 1),
            "best_low_coverage": float(np.max(u[low])),
        },
    }


def exact_audit(A: np.ndarray, name: str) -> None:
    n = len(A)
    reports = []
    for m in range((n + 1) // 2, n):
        q, selectors, records = ground_records(A, m)
        costs = np.array([r["cost"] for r in records], dtype=float)
        for threshold in (0, 4, 8):
            uniform = np.full(len(records), 1 / len(records))
            Z, Ccap, u = captured(records, threshold, uniform)
            covariance = float(np.mean(u * costs) - np.mean(u) * np.mean(costs))
            if m == n - 1:
                reports.append((threshold, Z, Ccap, covariance, sorted(Counter(zip(u, costs)).items())))

        for lam in (0.25, 0.5):
            # Budgets at the minimum and midpoint of the observed costs.
            for budget in sorted({float(costs.min()), float((costs.min() + costs.max()) / 2)}):
                result = coarea_audit(records, lam, budget)
                assert result["support"] <= 2

    print(name, {"n": n, "Q": q, "one_deletion_threshold_reports": reports})
    # Print representative nontrivial soft/coarea LP at one deletion.
    _, _, records = ground_records(A, n - 1)
    costs = np.array([r["cost"] for r in records], dtype=float)
    budget = float((costs.min() + costs.max()) / 2)
    print(name, "representative coarea", coarea_audit(records, 0.5, budget))


def abstract_gap_audit() -> None:
    # C=1.  A nearly-on-budget high-coverage point can rescue an
    # exponentially tiny low-cost point, whereas a fixed cost gap cannot.
    for L in (5.0, 10.0, 20.0):
        u_low = math.exp(-L)
        c_low = 0.0
        C = 1.0
        eps = u_low
        u_high = 1.0
        c_high = C + eps
        mu_low = (c_high - C) / (c_high - c_low)
        mu_high = 1 - mu_low
        Z_near = 1 / (mu_low / u_low + mu_high / u_high)
        assert Z_near > 0.49

        c_far = 2.0
        mu_low_far = (c_far - C) / (c_far - c_low)
        Z_far = 1 / (mu_low_far / u_low + (1 - mu_low_far) / u_high)
        assert Z_far <= 2.1 * u_low
        print(
            "abstract",
            {"L": L, "u_low": u_low, "near_gap": eps, "Z_near": Z_near, "Z_far": Z_far},
        )


def integer_budget_audit() -> None:
    """Audit the exact floor bound, including C<1 and nonintegral C."""
    for C in (0.0, 0.2, 1.0, 1.7, 7.999):
        lower_cost = 0
        upper_cost = math.floor(C) + 1
        if C == 0:
            low_mass = 1.0
        else:
            low_mass = 1 - C / upper_cost
        mu = np.asarray([low_mass, 1 - low_mass])
        costs = np.asarray([lower_cost, upper_cost], dtype=float)
        assert float(mu @ costs) <= C + 1e-12
        floor_budget = math.floor(C)
        exact_factor = (floor_budget + 1 - C) / (floor_budget + 1)
        assert math.isclose(float(mu[costs <= floor_budget].sum()), exact_factor)

        rounded_budget = math.ceil(C)
        rounded_mass = float(mu[costs <= rounded_budget].sum())
        assert rounded_mass + 1e-12 >= 1 / (rounded_budget + 1)
        print(
            "integer budget",
            {
                "C": C,
                "floor_factor": exact_factor,
                "rounded_budget": rounded_budget,
                "rounded_factor": 1 / (rounded_budget + 1),
            },
        )


def main() -> None:
    abstract_gap_audit()
    integer_budget_audit()
    for A, name in ((A6, "A6"), (A8, "A8"), (A9, "A9")):
        exact_audit(np.asarray(A, dtype=np.int64), name)
    print("PASS: switching/moment identities, two-ground LP, and loss-Laplace coarea recovery")


if __name__ == "__main__":
    main()
