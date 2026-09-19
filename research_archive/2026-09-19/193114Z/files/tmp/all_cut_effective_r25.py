#!/usr/bin/env python3
"""Exact finite audit of the arbitrary-cut effective-loss extension."""

from __future__ import annotations

import itertools
import math

import numpy as np
from scipy.optimize import linprog

from check_response_dual_r16 import A8, A9, spins
from verify_compatible_replacement_r12 import A6


def cut_records(A: np.ndarray, m: int):
    n = len(A)
    X = spins(n)
    parent_energies = np.einsum("bi,ij,bj->b", X, A, X)
    q = int(np.max(np.abs(parent_energies)))
    selectors = list(itertools.combinations(range(n), m))
    child_q = {}
    for S in selectors:
        B = A[np.ix_(S, S)]
        Y = spins(m)
        child_q[S] = int(np.max(np.abs(np.einsum("bi,ij,bj->b", Y, B, Y))))

    p_num = m * (m - 1)
    p_den = n * (n - 1)
    records = []
    for sigma in (-1, 1):
        for x, unsigned_parent in zip(X, parent_energies):
            parent = int(sigma * unsigned_parent)
            deficit = q - parent
            assert deficit >= 0
            rows = sigma * x * (x @ A)
            cost = int(rows @ rows)
            losses = []
            effective_numerators = []
            centered_numerators = []
            coefficient_good = []
            for S in selectors:
                xs = x[list(S)]
                child_score = int(sigma * xs @ A[np.ix_(S, S)] @ xs)
                loss = child_q[S] - child_score
                assert loss >= 0
                effective_num = p_den * loss - p_num * deficit
                centered_num = p_den * child_score - p_num * parent

                # Exact retained-deficit decomposition:
                # Q(A[S])-p_2 q = ell_eff + c_A(S,d)-p_2<A,d>.
                assert (
                    p_den * child_q[S] - p_num * q
                    == effective_num + centered_num
                )
                # Exact test of
                #   ell-p_2 Delta <= [(m/n)^(3/2)-p_2]q.
                # After adding p_2 q, square only when the rational left
                # side is positive.
                lhs_num = effective_num + p_num * q
                is_good = lhs_num <= 0 or (
                    lhs_num * lhs_num * n**3
                    <= q * q * m**3 * p_den * p_den
                )
                losses.append(loss)
                effective_numerators.append(effective_num)
                centered_numerators.append(centered_num)
                coefficient_good.append(is_good)

            records.append(
                {
                    "sigma": sigma,
                    "x": x.copy(),
                    "parent": parent,
                    "deficit": deficit,
                    "rows": rows.astype(int),
                    "cost": cost,
                    "losses": np.asarray(losses, dtype=int),
                    "effective_num": np.asarray(effective_numerators, dtype=int),
                    "centered_num": np.asarray(centered_numerators, dtype=int),
                    "coefficient_good": np.asarray(coefficient_good, dtype=bool),
                    "ground": deficit == 0,
                }
            )
    return q, selectors, p_num, p_den, records


def coverages(records, threshold_num: int = 0, coefficient_slack: bool = False):
    return np.asarray(
        [
            np.mean(record["coefficient_good"])
            if coefficient_slack
            else np.mean(record["effective_num"] <= threshold_num)
            for record in records
        ],
        dtype=float,
    )


def cost_lp(
    records, budget: float, ground_only: bool = False, coefficient_slack: bool = False
):
    indices = [i for i, r in enumerate(records) if r["ground"] or not ground_only]
    u = coverages(
        [records[i] for i in indices], coefficient_slack=coefficient_slack
    )
    c = np.asarray([records[i]["cost"] for i in indices], dtype=float)
    # Dummy action makes the probability simplex formulation always feasible.
    uu = np.r_[u, 0.0]
    cc = np.r_[c, 0.0]
    result = linprog(
        -uu,
        A_ub=np.asarray([uu * (cc - budget)]),
        b_ub=np.asarray([0.0]),
        A_eq=np.ones((1, len(uu))),
        b_eq=np.ones(1),
        bounds=(0, None),
        method="highs",
    )
    assert result.success
    prior = result.x[:-1]
    Z = float(prior @ u)
    captured_cost = float(prior @ (u * c) / Z) if Z else math.inf
    support = []
    for local_i in np.flatnonzero(prior > 1e-9):
        record = records[indices[int(local_i)]]
        support.append(
            {
                "prior": float(prior[local_i]),
                "captured": float(prior[local_i] * u[local_i] / Z),
                "u": float(u[local_i]),
                "cost": int(record["cost"]),
                "deficit": int(record["deficit"]),
            }
        )

    if Z:
        # Uniform exact extraction for a real budget: round upward.  The
        # captured mass on integer costs at most ceil(C) is >=1/(ceil(C)+1).
        integer_budget = int(math.ceil(budget))
        low = c <= integer_budget
        captured = prior * u / Z
        assert float(captured @ c) <= budget + 1e-8
        assert float(captured[low].sum()) + 1e-10 >= 1 / (integer_budget + 1)
        assert float(np.max(u[low])) + 1e-10 >= Z / (integer_budget + 1)
    return Z, captured_cost, support


def mgf_audit(A: np.ndarray, q: int, m: int, p_num: int, p_den: int, records) -> None:
    n = len(A)
    p = m / n
    p2 = p_num / p_den
    eps = p * p - p2
    prob_mode = math.comb(n, m) * p**m * (1 - p) ** (n - m)
    chi = -math.log(prob_mode)
    op = float(np.linalg.norm(A, ord=2))
    f_cross = float((n * n) // 4)
    lam = 0.1 / max(op, 1.0)
    assert 4 * lam * lam * op * op < 1
    for record in records:
        centered = record["centered_num"] / p_den
        a = float(np.max(lam * centered))
        # Stable finite log mean exp(lambda * centered).
        log_mgf = float(math.log(np.mean(np.exp(lam * centered - a))) + a)
        rhs = (
            chi
            + lam * eps * abs(record["parent"])
            + lam * lam * p * p * record["cost"]
            + lam * lam * f_cross / (1 - 4 * lam * lam * op * op)
        )
        assert log_mgf <= rhs + 1e-10


def audit(A: np.ndarray, name: str, budgets) -> None:
    A = np.asarray(A, dtype=np.int64)
    n = len(A)
    one_delete = None
    for m in range((n + 1) // 2, n):
        q, selectors, p_num, p_den, records = cut_records(A, m)
        mgf_audit(A, q, m, p_num, p_den, records)
        if m == n - 1:
            one_delete = records
    assert one_delete is not None

    for budget in budgets:
        all_value, all_cost, all_support = cost_lp(one_delete, budget, False)
        ground_value, ground_cost, ground_support = cost_lp(one_delete, budget, True)
        slack_value, slack_cost, slack_support = cost_lp(
            one_delete, budget, False, coefficient_slack=True
        )
        ground_slack_value, ground_slack_cost, ground_slack_support = cost_lp(
            one_delete, budget, True, coefficient_slack=True
        )
        print(
            name,
            {
                "budget": budget,
                "all_cut_Z": all_value,
                "all_cut_captured_cost": all_cost,
                "all_cut_support": all_support,
                "ground_Z": ground_value,
                "ground_captured_cost": ground_cost,
                "ground_support": ground_support,
                "coefficient_slack_all_cut_Z": slack_value,
                "coefficient_slack_all_cut_cost": slack_cost,
                "coefficient_slack_all_cut_support": slack_support,
                "coefficient_slack_ground_Z": ground_slack_value,
                "coefficient_slack_ground_cost": ground_slack_cost,
                "coefficient_slack_ground_support": ground_slack_support,
            },
        )


def main() -> None:
    audit(A6, "A6", (30,))
    audit(A8, "A8", (64,))
    audit(A9, "A9", (80, 88, 96, 104, 112))
    print("PASS: arbitrary-cut effective-loss identity, general mgf, LP, and extraction")


if __name__ == "__main__":
    main()
