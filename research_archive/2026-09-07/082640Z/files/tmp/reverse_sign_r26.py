#!/usr/bin/env python3
"""Finite audit for the Wave 26 centered Johnson sign lemma obstruction."""

from itertools import combinations

import numpy as np

from check_response_dual_r16 import A9, qnorm, spins


def selector_statistics(a: np.ndarray, m: int):
    n = len(a)
    q = qnorm(a)
    rho = m / n
    p2 = m * (m - 1) / (n * (n - 1))
    selectors = list(combinations(range(n), m))
    q_child = []
    for selected in selectors:
        block = a[np.ix_(selected, selected)]
        q_child.append(qnorm(block))

    records = []
    for x in spins(n):
        parent_unsigned = int(x @ a @ x)
        row_square = int(np.sum((a @ x) ** 2))
        for sigma in (-1, 1):
            parent = sigma * parent_unsigned
            x_values = []
            y_values = []
            for selected, child_q in zip(selectors, q_child):
                xs = x[list(selected)]
                child = int(sigma * xs @ a[np.ix_(selected, selected)] @ xs)
                x_values.append(child - p2 * parent)
                y_values.append(child_q - rho**1.5 * q)
            xv = np.asarray(x_values)
            yv = np.asarray(y_values)
            assert abs(float(np.mean(xv))) < 1e-12
            records.append(
                {
                    "row_square": row_square,
                    "positive": int(np.sum(xv > 1e-12)),
                    "negative": int(np.sum(xv < -1e-12)),
                    "coverage": int(np.sum(yv - xv <= 1e-12)),
                    "max_x": float(np.max(xv)),
                    "min_y": float(np.min(yv)),
                    "max_y": float(np.max(yv)),
                }
            )
    return selectors, records


def main():
    for m, expected_positive in ((8, 5), (7, 22)):
        selectors, records = selector_statistics(A9.astype(np.int64), m)
        minimum_row = min(record["row_square"] for record in records)
        witnesses = [
            record
            for record in records
            if record["row_square"] == minimum_row
            and record["positive"] == expected_positive
            and record["coverage"] == 0
        ]
        assert minimum_row == 16
        assert witnesses
        witness = witnesses[0]
        print(
            {
                "m": m,
                "selectors": len(selectors),
                "minimum_R2": minimum_row,
                "positive_X": witness["positive"],
                "negative_X": witness["negative"],
                "effective_zero_coverage": witness["coverage"],
                "max_X": witness["max_x"],
                "Y_range": (witness["min_y"], witness["max_y"]),
            }
        )
    print("centered Johnson finite obstruction: all checks passed")


if __name__ == "__main__":
    main()
