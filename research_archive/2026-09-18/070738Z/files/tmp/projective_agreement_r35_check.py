#!/usr/bin/env python3
"""Exact finite checks for the Wave 35 anchored agreement decoder.

This is only an audit of the elementary decoder and small A6/A8/A9 data.
It is not evidence for an asymptotic agreement or replacement theorem.
"""

from fractions import Fraction
from itertools import combinations, product
from math import floor

import numpy as np

from coherent_signatures_r29_search import A6, A8, A9, x_of


def child_labels(a: np.ndarray, selector: tuple[int, ...], allowance: int):
    idx = np.array(selector)
    child = a[np.ix_(idx, idx)]
    labels = [np.array((1,) + tail, dtype=np.int64)
              for tail in product((-1, 1), repeat=len(selector) - 1)]
    values = [abs(int(y @ child @ y)) for y in labels]
    optimum = max(values)
    return [y for y, value in zip(labels, values)
            if optimum - value <= allowance]


def anchored_decoder(selectors, labels, n, anchor=(0,)):
    """Return a plurality word, its mean error, and degree-corrected conflict."""
    count = len(selectors)
    assert count
    z = np.ones(n, dtype=np.int64)
    conflict = Fraction(0)
    anchor = set(anchor)
    for i in range(n):
        if i in anchor:
            assert all(labels[s][s.index(i)] == 1 for s in selectors)
            continue
        plus = sum(i in s and labels[s][s.index(i)] == 1 for s in selectors)
        minus = sum(i in s and labels[s][s.index(i)] == -1 for s in selectors)
        degree = plus + minus
        if not degree:
            continue
        z[i] = 1 if plus >= minus else -1
        # 2 a_i b_i / p_i after cancelling the common |G| denominator.
        conflict += Fraction(2 * plus * minus, count * degree)
    total_error = sum(
        int(np.count_nonzero(labels[s] != z[list(s)])) for s in selectors
    )
    mean_error = Fraction(total_error, count)
    assert mean_error <= conflict
    return z, mean_error, conflict


def audit_matrix(a: np.ndarray, name: str, m: int):
    n = len(a)
    full = [x_of(mask, n) for mask in range(1 << (n - 1))]
    q = max(abs(int(x @ a @ x)) for x in full)
    allowance = floor(((m / n) ** 1.5 - (m / n) ** 2) * q + 1e-12)
    all_selectors = list(combinations(range(n), m))
    lists = {s: child_labels(a, s, allowance) for s in all_selectors}

    degrees = []
    rows = []
    for x in full:
        degree = 0
        for s in all_selectors:
            target = x[list(s)]
            if any(min(np.count_nonzero(y != target),
                       np.count_nonzero(y != -target)) == 0 for y in lists[s]):
                degree += 1
        degrees.append(degree)
        rows.append(int((a @ x) @ (a @ x)))

    best_index = max(range(len(full)), key=lambda j: degrees[j])
    best = full[best_index]
    anchored = [s for s in all_selectors if 0 in s]
    chosen = {}
    for s in anchored:
        target = best[list(s)]
        # Since global coordinate 0 is first in s and is +1, these are the
        # canonical anchor-oriented representatives.
        chosen[s] = min(lists[s], key=lambda y: int(np.count_nonzero(y != target)))
    median, mean_error, conflict = anchored_decoder(anchored, chosen, n)

    cap = 2 * n * (n - 1)
    assert max(rows) <= cap
    assert int((a @ median) @ (a @ median)) <= cap
    return {
        "name": name,
        "m": m,
        "allowance": allowance,
        "max_degree": Fraction(max(degrees), len(all_selectors)),
        "row_at_chosen_max": rows[best_index],
        "median_row": int((a @ median) @ (a @ median)),
        "mean_error": mean_error,
        "conflict": conflict,
        "row_range": (min(rows), max(rows)),
    }


def synthetic_zero_degree_audit():
    """Exercise the p_i=0 convention in the decoder."""
    selectors = [(0, 1), (0, 1)]
    labels = {
        selectors[0]: np.array([1, 1], dtype=np.int64),
    }
    # Dict keys coincide, so use the same label; coordinates 2,3 never occur.
    z, mean_error, conflict = anchored_decoder(selectors, labels, 4)
    assert tuple(z) == (1, 1, 1, 1)
    assert mean_error == conflict == 0


if __name__ == "__main__":
    synthetic_zero_degree_audit()
    rows = [
        audit_matrix(A6, "A6", 3),
        audit_matrix(A8, "A8", 5),
        audit_matrix(A9, "A9", 6),
    ]
    for row in rows:
        print(row)
    print("PASS: anchored decoder inequality and finite center tables")
