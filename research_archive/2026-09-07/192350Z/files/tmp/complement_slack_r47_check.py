#!/usr/bin/env python3
"""Exact finite checks for Wave 47A complement-slack replacement walls."""

from __future__ import annotations

import itertools
import sys
from collections import Counter
from fractions import Fraction

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from envelope_block_cover_r27 import A6, A8, A9


def projective_spins(n: int):
    for tail in itertools.product((-1, 1), repeat=n - 1):
        yield np.asarray((1,) + tail, dtype=np.int64)


def qvalue(a: np.ndarray) -> int:
    return max(abs(int(x @ a @ x)) for x in projective_spins(len(a)))


def switch(a: np.ndarray, vertices: tuple[int, ...]) -> np.ndarray:
    eps = np.ones(len(a), dtype=np.int64)
    eps[list(vertices)] = -1
    return eps[:, None] * a * eps[None, :]


def flip_edges(a: np.ndarray, edges: tuple[tuple[int, int], ...]) -> np.ndarray:
    b = a.copy()
    for i, j in edges:
        b[i, j] *= -1
        b[j, i] *= -1
    return b


def complement_signing(a: np.ndarray, s: tuple[int, ...]) -> np.ndarray:
    b = -a.copy()
    idx = np.asarray(s)
    b[np.ix_(idx, idx)] = a[np.ix_(idx, idx)]
    return b


def build_states(a: np.ndarray, m: int):
    n = len(a)
    q = qvalue(a)
    selectors = list(itertools.combinations(range(n), m))
    states = []
    for x in projective_spins(n):
        row = int((a @ x) @ (a @ x))
        raw = int(x @ a @ x)
        for sigma in (-1, 1):
            scores = []
            for s in selectors:
                idx = np.asarray(s)
                child = sigma * int(x[idx] @ a[np.ix_(idx, idx)] @ x[idx])
                scores.append(2 * child - sigma * raw)
            alpha = Fraction(sum(v >= q for v in scores), len(selectors))
            if alpha:
                states.append((row, alpha, sigma, x.copy(), tuple(scores)))
    return q, selectors, states


def pareto_representatives(states):
    pareto = [
        d
        for d in states
        if not any(
            e[0] <= d[0]
            and e[1] >= d[1]
            and (e[0] < d[0] or e[1] > d[1])
            for e in states
        )
    ]
    reps = {}
    for d in pareto:
        row, alpha, _sigma, _x, scores = d
        q_dummy = None
        # The score multiset itself distinguishes all audited Pareto types once
        # the active threshold is supplied by audit().
        reps.setdefault((row, alpha, scores), d)
    return pareto, reps


def audit(a: np.ndarray, m: int, name: str):
    n = len(a)
    q, selectors, states = build_states(a, m)
    pareto = [
        d
        for d in states
        if not any(
            e[0] <= d[0]
            and e[1] >= d[1]
            and (e[0] < d[0] or e[1] > d[1])
            for e in states
        )
    ]
    reps = {}
    for d in pareto:
        row, alpha, _sigma, _x, scores = d
        slacks = tuple(sorted(v - q for v in scores if v >= q))
        reps.setdefault((row, alpha, slacks), d)

    incidence_profiles = Counter()
    for d in pareto:
        _row, _alpha, _sigma, _x, scores = d
        for s, score in zip(selectors, scores):
            if score < q:
                continue
            t = tuple(sorted(set(range(n)) - set(s)))
            internal_t = tuple(itertools.combinations(t, 2))
            b_s = complement_signing(a, s)
            small_flip = flip_edges(a, internal_t)
            assert np.array_equal(b_s, switch(small_flip, t))
            cap = qvalue(b_s)
            gamma = score - q
            assert 0 <= gamma <= cap - q <= 2 * len(t) * (len(t) - 1)
            incidence_profiles[(gamma, cap - q, cap - score)] += 1

    combination_profiles = []
    combination_signatures_by_type = {}
    for d in pareto:
        row, alpha, sigma, x, scores = d
        kind = (row, alpha, tuple(sorted(v - q for v in scores if v >= q)))
        active = []
        for s, score in zip(selectors, scores):
            if score < q:
                continue
            t = tuple(sorted(set(range(n)) - set(s)))
            internal_t = tuple(itertools.combinations(t, 2))
            b_s = complement_signing(a, s)
            small_flip = flip_edges(a, internal_t)
            assert np.array_equal(b_s, switch(small_flip, t))
            cap = qvalue(b_s)
            gamma = score - q
            assert 0 <= gamma <= cap - q <= 2 * len(t) * (len(t) - 1)
            y = x.copy()
            y[list(t)] *= -1
            assert sigma * int(y @ small_flip @ y) == score
            active.append((t, internal_t, y, score))

        # At m=n-2, combine the natural missing-pair flips attached to every
        # active incidence of this one Pareto column.  This is the most literal
        # common replacement suggested by the active family.
        active_edges = tuple(sorted({e for _t, es, _y, _score in active for e in es}))
        combined = flip_edges(a, active_edges)
        selected_energies = tuple(
            sigma * int(y @ combined @ y) for _t, _es, y, _score in active
        )
        cap = qvalue(combined)
        shield_count = sum(
            sigma2 * int(z @ combined @ z) >= q
            for z in projective_spins(n)
            for sigma2 in (-1, 1)
        )

        subset_caps = []
        for mask in range(1 << len(active_edges)):
            f = tuple(e for j, e in enumerate(active_edges) if (mask >> j) & 1)
            subset_caps.append(qvalue(flip_edges(a, f)))
        assert min(subset_caps) == q
        signature = (
            len(active_edges),
            selected_energies,
            cap,
            shield_count,
            tuple(sorted(Counter(subset_caps).items())),
        )
        combination_signatures_by_type.setdefault(kind, Counter())[signature] += 1

    for kind, signatures in combination_signatures_by_type.items():
        for signature, multiplicity in signatures.items():
            edge_count, selected_energies, cap, shield_count, subset_hist = signature
            combination_profiles.append(
                {
                    "pareto_type": str(kind),
                    "multiplicity": multiplicity,
                    "active_missing_edge_count": edge_count,
                    "combined_selected_energies": selected_energies,
                    "combined_cap": cap,
                    "combined_shield_states_at_least_q": shield_count,
                    "subset_cap_histogram": dict(subset_hist),
                }
            )

    # Universal star-switch identity at m=n-1: every active score is exactly
    # q, independently of scalar price or Pareto optimality.
    q_star, star_selectors, star_states = build_states(a, n - 1)
    assert q_star == q
    star_slacks = []
    for _row, _alpha, sigma, x, scores in star_states:
        for s, score in zip(star_selectors, scores):
            if score >= q:
                t = tuple(sorted(set(range(n)) - set(s)))
                assert len(t) == 1
                b_s = complement_signing(a, s)
                assert np.array_equal(b_s, switch(a, t))
                assert qvalue(b_s) == q
                star_slacks.append(score - q)
    assert star_slacks and set(star_slacks) == {0}

    print(
        {
            "name": name,
            "n": n,
            "m": m,
            "q": q,
            "incidence_(Gamma,cap_excess,unseen_excess)_counts": dict(
                sorted(incidence_profiles.items())
            ),
            "pareto_common_replacement_profiles": combination_profiles,
            "m=n-1_active_incidence_count": len(star_slacks),
            "m=n-1_active_slacks": sorted(set(star_slacks)),
        }
    )


def main():
    audit(A6, 5, "A6")
    audit(A8, 6, "A8")
    audit(A9, 7, "A9")
    print("PASS complement_slack_r47_check")


if __name__ == "__main__":
    main()
