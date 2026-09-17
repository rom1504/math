#!/usr/bin/env python3
"""Exact finite checks for the structured ground-flip incidence in Wave 38."""

from __future__ import annotations

import itertools
import math
import sys

import numpy as np
from scipy.optimize import linprog

sys.path.insert(0, "/home/math/quadra/tmp")
from envelope_block_cover_r27 import A6, A8, A9


def projective_spins(n: int):
    for tail in itertools.product((-1, 1), repeat=n - 1):
        yield np.asarray((1,) + tail, dtype=np.int64)


def q_value(a: np.ndarray) -> int:
    return max(abs(int(x @ a @ x)) for x in projective_spins(len(a)))


def fractional_cover(inc: np.ndarray) -> float | None:
    if np.any(inc.sum(axis=1) == 0):
        return None
    sol = linprog(
        np.ones(inc.shape[1]),
        A_ub=-inc,
        b_ub=-np.ones(inc.shape[0]),
        bounds=(0.0, None),
        method="highs",
    )
    assert sol.success
    return float(sol.fun)


def flipped_complement(a: np.ndarray, selector: tuple[int, ...]) -> np.ndarray:
    n = len(a)
    inside = np.zeros(n, dtype=bool)
    inside[list(selector)] = True
    keep = np.outer(inside, inside)
    b = np.where(keep, a, -a)
    np.fill_diagonal(b, 0)
    return b


def audit(a: np.ndarray, m: int, name: str) -> None:
    n = len(a)
    q = q_value(a)
    selectors = list(itertools.combinations(range(n), m))
    states = []
    for x in projective_spins(n):
        for sigma in (-1, 1):
            signed_degrees = sigma * x * (a @ x)
            energy = int(signed_degrees.sum())
            row = int(signed_degrees @ signed_degrees)
            positive = tuple(np.flatnonzero(signed_degrees > 0).tolist())
            states.append((x, sigma, signed_degrees, energy, row, positive))

    ground = np.zeros((len(selectors), len(states)), dtype=np.int8)
    incidence = np.zeros_like(ground)
    flip_caps = []

    for si, selector in enumerate(selectors):
        b = flipped_complement(a, selector)
        qb = q_value(b)
        flip_caps.append(qb)
        selector_set = set(selector)
        for di, (x, sigma, degrees, energy, _, positive) in enumerate(states):
            flip_energy = sigma * int(x @ b @ x)
            incidence[si, di] = flip_energy >= q
            ground[si, di] = flip_energy == qb
            if ground[si, di]:
                assert incidence[si, di]
                assert set(positive).issubset(selector_set)

                # Full ground compatibility is equivalent to nonnegative
                # switched cut weight for every vertex subset.
                signed_edges = sigma * np.outer(x, x) * a
                switched_edges = sigma * np.outer(x, x) * b
                for bits in itertools.product((0, 1), repeat=n):
                    u = np.asarray(bits, dtype=bool)
                    cut = np.logical_xor.outer(u, u)
                    switched_cut = int(switched_edges[cut].sum() // 2)
                    parent_cut = int(signed_edges[cut].sum() // 2)
                    internal_cut = int(
                        signed_edges[
                            cut & np.outer(
                                np.isin(np.arange(n), selector),
                                np.isin(np.arange(n), selector),
                            )
                        ].sum()
                        // 2
                    )
                    assert switched_cut >= 0
                    assert switched_cut == 2 * internal_cut - parent_cut

    assert np.all(ground.sum(axis=1) >= 1)
    assert np.all(ground <= incidence)

    ground_degrees = ground.sum(axis=0)
    for di, (*_, positive) in enumerate(states):
        k = len(positive)
        upper = math.comb(n - k, m - k) if k <= m else 0
        assert int(ground_degrees[di]) <= upper

    row = np.asarray([state[4] for state in states])
    min_ground_row = np.min(np.where(ground, row[None, :], 10**18), axis=1)
    min_inc_row = np.min(np.where(incidence, row[None, :], 10**18), axis=1)
    mean_cap = n * (n - 1)
    projective_caps = sorted(set(row.tolist()))
    cover_curve = []
    for cap in projective_caps:
        cols = row <= cap
        if np.all(ground[:, cols].sum(axis=1) > 0):
            cover_curve.append((int(cap), fractional_cover(ground[:, cols])))
            break

    active = np.flatnonzero(ground_degrees)
    by_degree = sorted(
        (
            int(ground_degrees[di]),
            int(states[di][4]),
            len(states[di][5]),
            int(states[di][3]),
        )
        for di in active
    )
    max_ground_degree = max(v[0] for v in by_degree)
    max_profiles = [v for v in by_degree if v[0] == max_ground_degree]

    print(
        name,
        {
            "n": n,
            "m": m,
            "q": q,
            "flip_cap_hist": {
                value: flip_caps.count(value) for value in sorted(set(flip_caps))
            },
            "ground_witness_count_range": (
                int(ground.sum(axis=1).min()),
                int(ground.sum(axis=1).max()),
            ),
            "min_ground_row_hist": {
                int(value): int(np.sum(min_ground_row == value))
                for value in sorted(set(min_ground_row.tolist()))
            },
            "min_incidence_row_hist": {
                int(value): int(np.sum(min_inc_row == value))
                for value in sorted(set(min_inc_row.tolist()))
            },
            "ground_mean_cap_misses": int(np.sum(min_ground_row > mean_cap)),
            "first_full_ground_cap_and_numeric_tau": cover_curve[0],
            "unrestricted_ground_tau_numeric": fractional_cover(ground),
            "max_state_ground_degree": max_ground_degree,
            "max_degree_profiles_degree_row_pos_energy": max_profiles[:8],
        },
    )


def main() -> None:
    audit(A6, 5, "A6")
    audit(A8, 6, "A8")
    audit(A9, 7, "A9")
    print("PASS ground_flip_incidence_r38_check")


if __name__ == "__main__":
    main()
