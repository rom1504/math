#!/usr/bin/env python3
"""Exact finite audit of row-priced complement-flip witnesses (Wave 39)."""

from __future__ import annotations

import itertools
from fractions import Fraction

import numpy as np
from scipy.optimize import linprog

from envelope_block_cover_r27 import A6, A8, A9


def projective_spins(n: int):
    for tail in itertools.product((-1, 1), repeat=n - 1):
        yield np.asarray((1,) + tail, dtype=np.int64)


def q_value(a: np.ndarray) -> int:
    return max(abs(int(x @ a @ x)) for x in projective_spins(len(a)))


def flipped_complement(a: np.ndarray, selector: tuple[int, ...]) -> np.ndarray:
    inside = np.zeros(len(a), dtype=bool)
    inside[list(selector)] = True
    b = np.where(np.outer(inside, inside), a, -a)
    np.fill_diagonal(b, 0)
    return b


def critical_price(profiles: list[tuple[int, int]], q: int) -> Fraction:
    """Least lambda for which a high-energy profile is row-price optimal."""
    high = [(ell, row) for ell, row in profiles if ell >= q]
    low = [(ell, row) for ell, row in profiles if ell < q]
    answer: Fraction | None = None
    for ell_h, row_h in high:
        needed = Fraction(0)
        possible = True
        for ell_l, row_l in low:
            if ell_h <= ell_l:
                if row_h > row_l:
                    possible = False
                    break
                continue
            needed = max(needed, Fraction(row_h - row_l, ell_h - ell_l))
        if possible and (answer is None or needed < answer):
            answer = needed
    assert answer is not None
    return answer


def audit(a: np.ndarray, m: int, name: str) -> None:
    n = len(a)
    q = q_value(a)
    spins = list(projective_spins(n))
    rows = [int(x @ a @ a @ x) for x in spins]
    selectors = list(itertools.combinations(range(n), m))
    prices = []
    chosen_rows = []
    chosen_energies = []
    for selector in selectors:
        b = flipped_complement(a, selector)
        profiles = sorted(
            set((abs(int(x @ b @ x)), row) for x, row in zip(spins, rows))
        )
        price = critical_price(profiles, q)
        prices.append(price)
        vals = [(price * ell - row, ell, row) for ell, row in profiles]
        best = max(v[0] for v in vals)
        eligible = [(ell, row) for value, ell, row in vals if value == best and ell >= q]
        assert eligible
        chosen_energies.append(max(ell for ell, _ in eligible))
        chosen_rows.append(min(row for _, row in eligible))

        # The global minimizer of R-lambda L obeys R <= n(n-1)+lambda L.
        for ell, row in eligible:
            assert row <= n * (n - 1) + price * ell

    def hist(values):
        return {str(v): values.count(v) for v in sorted(set(values))}

    common_price = max(prices)
    oriented = [(x, sigma) for x in spins for sigma in (-1, 1)]
    optimal_incidence = np.zeros((len(selectors), len(oriented)), dtype=np.int8)
    for si, selector in enumerate(selectors):
        b = flipped_complement(a, selector)
        scores = [
            common_price * sigma * int(x @ b @ x) - int(x @ a @ a @ x)
            for x, sigma in oriented
        ]
        best = max(scores)
        for di, ((x, sigma), score) in enumerate(zip(oriented, scores)):
            ell = sigma * int(x @ b @ x)
            if score == best and ell >= q:
                optimal_incidence[si, di] = 1
    assert np.all(optimal_incidence.sum(axis=1) > 0)
    active = np.flatnonzero(optimal_incidence.sum(axis=0))
    lp = linprog(
        np.ones(len(active)),
        A_ub=-optimal_incidence[:, active],
        b_ub=-np.ones(len(selectors)),
        bounds=(0.0, None),
        method="highs",
    )
    assert lp.success

    print(
        name,
        {
            "q": q,
            "critical_price_hist": hist(prices),
            "max_critical_price": str(max(prices)),
            "critical_energy_hist": hist(chosen_energies),
            "critical_row_hist": hist(chosen_rows),
            "common_price_active_oriented_states": len(active),
            "common_price_fractional_cover_numeric": float(lp.fun),
        },
    )


def main() -> None:
    audit(A6, 5, "A6")
    audit(A8, 6, "A8")
    audit(A9, 7, "A9")
    print("PASS complement_row_price_r39")


if __name__ == "__main__":
    main()
