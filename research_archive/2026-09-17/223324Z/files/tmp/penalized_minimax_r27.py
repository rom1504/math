#!/usr/bin/env python3
"""Finite audit of the Wave 27 penalized all-cut mean dual."""

from itertools import combinations

import numpy as np

from check_response_dual_r16 import A8, A9, qnorm, spins
from verify_compatible_replacement_r12 import A6


def audit(a: np.ndarray, m: int) -> None:
    a = a.astype(np.int64)
    n = len(a)
    q = qnorm(a)
    rho = m / n
    p2 = m * (m - 1) / (n * (n - 1))
    selectors = list(combinations(range(n), m))
    child_q = np.asarray(
        [qnorm(a[np.ix_(selected, selected)]) for selected in selectors],
        dtype=float,
    )

    # A deterministic nonuniform selector law, plus the uniform law.
    raw = np.arange(1, len(selectors) + 1, dtype=float) ** 2
    laws = (np.full(len(selectors), 1 / len(selectors)), raw / raw.sum())
    x_reps = spins(n)

    for weights in laws:
        inclusion = np.zeros((n, n), dtype=float)
        for mass, selected in zip(weights, selectors):
            inclusion[np.ix_(selected, selected)] += mass
        centered = a * (inclusion - p2)
        np.fill_diagonal(centered, 0.0)
        y_bar = float(weights @ child_q - rho**1.5 * q)

        records = []
        for x in x_reps:
            parent = int(x @ a @ x)
            row_square = int(np.sum((a @ x) ** 2))
            x_unsigned = float(x @ centered @ x)
            for sigma in (-1, 1):
                effective = []
                for selected, q_s in zip(selectors, child_q):
                    xs = x[list(selected)]
                    child = int(sigma * xs @ a[np.ix_(selected, selected)] @ xs)
                    x_value = child - p2 * sigma * parent
                    effective.append(q_s - rho**1.5 * q - x_value)
                records.append(
                    (row_square, float(weights @ effective), sigma * x_unsigned)
                )

        costs = sorted({record[0] for record in records})
        for budget in (costs[0], costs[len(costs) // 2], costs[-1]):
            eligible = [record for record in records if record[0] <= budget]
            direct = min(record[1] for record in eligible)
            dual_expression = y_bar - max(record[2] for record in eligible)
            assert abs(direct - dual_expression) < 1e-10

        for price in (0.0, 0.01, 0.1):
            budget = costs[len(costs) // 2]
            direct = min(
                mean_loss + price * (cost - budget)
                for cost, mean_loss, _ in records
            )
            formula = y_bar - price * budget - max(
                x_value - price * cost for cost, _, x_value in records
            )
            assert abs(direct - formula) < 1e-10

        if np.allclose(weights, laws[0]):
            assert np.max(np.abs(centered)) < 1e-12
            means = np.asarray([record[1] for record in records])
            assert np.max(np.abs(means - y_bar)) < 1e-10


if __name__ == "__main__":
    for matrix in (A6, A8, A9):
        audit(matrix, len(matrix) - 1)
        audit(matrix, len(matrix) - 2)
    print("penalized all-cut mean dual: all finite checks passed")
