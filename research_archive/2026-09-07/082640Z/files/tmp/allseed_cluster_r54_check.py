#!/usr/bin/env python3
"""Reproducible checks for the Wave 54 all-seed/cluster audit."""

from __future__ import annotations

import itertools
import math
import sys

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from box_triple_r50_check import A  # exact order-ten minimizer


def entropy(x: float) -> float:
    if x <= 0.0 or x >= 1.0:
        return 0.0
    return -x * math.log(x) - (1.0 - x) * math.log(1.0 - x)


def ball_rate(p: float, beta: float) -> float:
    """Exponential rate of the constant-weight ball |S cap T| >= beta*n."""
    if beta <= p * p:
        return entropy(p)
    delta = p - beta
    return p * entropy(delta / p) + (1.0 - p) * entropy(delta / (1.0 - p))


def required_rate(p: float, alpha: float, beta: float) -> float:
    """Exponential rate of the zero-epsilon RHS in (10.1294)."""
    return (
        p * entropy(alpha / p)
        + (1.0 - alpha) * entropy((p - alpha) / (1.0 - alpha))
        - beta * entropy(alpha / beta)
    )


def entropy_audit() -> None:
    min_gap = float("inf")
    tested = 0
    # A deterministic grid checks both sides of beta=p^2.  The proof is in
    # the accompanying memo; this is only a floating-point regression check.
    for ip in range(51, 100):
        p = ip / 100.0
        for ib in range(1, 100):
            beta = p * ib / 100.0
            for ia in range(1, ib):
                alpha = p * ia / 100.0
                gap = required_rate(p, alpha, beta) - ball_rate(p, beta)
                assert gap >= -2.0e-12
                min_gap = min(min_gap, gap)
                tested += 1

    examples = []
    for p, beta in ((0.6, 0.45), (0.7, 0.6), (0.85, 0.8)):
        assert beta > p * p
        alpha_star = (beta - p * p) / (1.0 + beta - 2.0 * p)
        gap = required_rate(p, alpha_star, beta) - ball_rate(p, beta)
        assert abs(gap) < 2.0e-12
        examples.append((p, beta, alpha_star, gap))

    p, alpha, beta = 0.6, 0.3, 0.45
    e = required_rate(p, alpha, beta)
    v = ball_rate(p, beta)
    assert e > v
    print("entropy_grid", tested, "minimum_gap", f"{min_gap:.12g}")
    for row in examples:
        print("equality_case", *(f"{x:.12g}" for x in row))
    print("ledger_example_rates", f"E={e:.12g}", f"V={v:.12g}", f"gap={e-v:.12g}")


def exact_ledger_example() -> None:
    """Compare (10.1294)'s distinct-partner demand with the whole ball."""
    first_impossible = None
    rows = []
    for n in range(20, 1001, 20):
        m = 3 * n // 5
        ell = 3 * n // 10
        s = 9 * n // 20
        b_ell = math.comb(m - 2, ell - 2) * math.comb(n - ell - 2, m - ell)
        diagonal = math.comb(m, ell)
        # This is the epsilon=0 RHS.  Any positive epsilon only increases it.
        demand_num = b_ell - diagonal
        demand_den = math.comb(s, ell)
        # Distinct partners only: omit j=m, which is the self-loop T=S.
        whole_ball = sum(
            math.comb(m, j) * math.comb(n - m, m - j)
            for j in range(s, m)
        )
        impossible = demand_num > whole_ball * demand_den
        if impossible and first_impossible is None:
            first_impossible = n
        if n in (20, 100, 160, 200, 400, 1000):
            log_ratio_per_n = (
                math.log(demand_num) - math.log(demand_den) - math.log(whole_ball)
            ) / n
            rows.append((n, impossible, log_ratio_per_n))
    assert first_impossible == 160
    print("ledger_example_first_exact_impossible_n", first_impossible)
    for n, impossible, rate in rows:
        print("ledger_exact", n, "impossible", impossible, "log_ratio_per_n", f"{rate:.12g}")


def projective_spins(n: int):
    for tail in itertools.product((-1, 1), repeat=n - 1):
        yield np.asarray((1,) + tail, dtype=np.int64)


def qnorm(matrix: np.ndarray) -> int:
    return max(abs(int(x @ matrix @ x)) for x in projective_spins(len(matrix)))


def arbitrary_seed_star_a10() -> None:
    """Exhaust the arbitrary-seed H=8(m-1) theorem on A10 at m=6."""
    n = len(A)
    m = 6
    k = n - m
    selectors = list(itertools.combinations(range(n), m))
    caps = {s: qnorm(A[np.ix_(s, s)]) for s in selectors}

    max_cap_lipschitz = 0
    max_budget_ratio = 0.0
    max_formula_error = 0
    minimum_assigned_star = m * k
    incidences = 0

    for selector in selectors:
        child = A[np.ix_(selector, selector)]
        outside = tuple(i for i in range(n) if i not in selector)
        local_index = {v: i for i, v in enumerate(selector)}
        grounds = []
        for y in projective_spins(m):
            value = int(y @ child @ y)
            if abs(value) == caps[selector]:
                grounds.append((y, value))

        for y, value in grounds:
            sigma = 1 if value > 0 else -1
            for tail in itertools.product((-1, 1), repeat=k):
                z = np.empty(n, dtype=np.int64)
                z[list(selector)] = y
                z[list(outside)] = tail
                good_by_center: dict[tuple[int, ...], set[tuple[int, ...]]] = {}
                deficit_sum = 0

                for i in selector:
                    ii = local_index[i]
                    core = tuple(u for u in selector if u != i)
                    r_i = sigma * int(y[ii]) * sum(
                        int(A[i, u]) * int(y[local_index[u]]) for u in core
                    )
                    assert 0 <= r_i <= caps[selector] / 2

                    for j in outside:
                        exchanged = tuple(sorted(core + (j,)))
                        cap_jump = abs(caps[exchanged] - caps[selector])
                        max_cap_lipschitz = max(max_cap_lipschitz, cap_jump)
                        assert cap_jump <= 2 * (m - 1)

                        best_deficit = None
                        best_center = None
                        for flip in (False, True):
                            center = z.copy()
                            if flip:
                                center[j] *= -1
                            x = center[list(exchanged)]
                            candidate = abs(int(x @ A[np.ix_(exchanged, exchanged)] @ x))
                            deficit = caps[exchanged] - candidate
                            if best_deficit is None or deficit < best_deficit:
                                best_deficit = deficit
                                best_center = tuple(map(int, center))

                        assert best_deficit is not None and best_center is not None
                        h_ji = sum(
                            int(A[j, u]) * int(y[local_index[u]]) for u in core
                        )
                        formula = (
                            caps[exchanged]
                            - caps[selector]
                            + 2 * (r_i - abs(h_ji))
                        )
                        max_formula_error = max(max_formula_error, abs(best_deficit - formula))
                        assert best_deficit == formula and best_deficit >= 0
                        deficit_sum += best_deficit

                        if best_deficit <= 8 * (m - 1):
                            good_by_center.setdefault(best_center, set()).add(exchanged)

                budget = 4 * m * k * (m - 1)
                assert deficit_sum <= budget
                max_budget_ratio = max(max_budget_ratio, deficit_sum / budget)
                assigned_star = max(map(len, good_by_center.values()))
                assert assigned_star >= m * k / (2 * (k + 1))
                minimum_assigned_star = min(minimum_assigned_star, assigned_star)
                incidences += 1

    print(
        "a10_arbitrary_seed_star",
        "incidences", incidences,
        "min_assigned_star", minimum_assigned_star,
        "max_cap_lipschitz", max_cap_lipschitz,
        "max_budget_ratio", f"{max_budget_ratio:.12g}",
        "max_formula_error", max_formula_error,
    )


if __name__ == "__main__":
    entropy_audit()
    exact_ledger_example()
    arbitrary_seed_star_a10()
