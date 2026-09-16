#!/usr/bin/env python3
"""Exact checks for Wave 46 noncentral complement-slack switching."""

from __future__ import annotations

import itertools
import math
import sys
from fractions import Fraction

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from anchored_conflict_r40_check import A5
from envelope_block_cover_r27 import A6, A8, A9


def projective_spins(n: int):
    for tail in itertools.product((-1, 1), repeat=n - 1):
        yield np.asarray((1,) + tail, dtype=np.int64)


def qvalue(a: np.ndarray) -> int:
    return max(abs(int(x @ a @ x)) for x in projective_spins(len(a)))


def canonical(x: np.ndarray) -> tuple[int, ...]:
    y = x.copy()
    if y[0] < 0:
        y *= -1
    return tuple(int(v) for v in y)


def eps_moments(n: int, k: int) -> tuple[Fraction, Fraction]:
    # Exact averages of eps_i eps_j and eps_i eps_j eps_l eps_r on a k-layer.
    theta = Fraction(n * (n - 1) - 4 * k * (n - k), n * (n - 1))
    mu4 = sum(
        Fraction((-1) ** j * math.comb(k, j) * math.comb(n - k, 4 - j),
                 math.comb(n, 4))
        for j in range(5)
        if j <= k and 4 - j <= n - k
    )
    return theta, mu4


def exact_variance_formula(
    b: np.ndarray, x: np.ndarray, sigma: int, k: int
) -> tuple[Fraction, Fraction, Fraction, int, int]:
    n = len(b)
    theta, mu4 = eps_moments(n, k)
    l0 = sigma * int(x @ b @ x)
    rb = int((b @ x) @ (b @ x))
    predicted = (
        (mu4 - theta * theta) * l0 * l0
        + 4 * (theta - mu4) * rb
        + 2 * (1 - 2 * theta + mu4) * n * (n - 1)
    )
    values = []
    for u in itertools.combinations(range(n), k):
        y = x.copy()
        y[list(u)] *= -1
        values.append(sigma * int(y @ b @ y))
    mean = Fraction(sum(values), len(values))
    empirical = sum((Fraction(v) - mean) ** 2 for v in values) / len(values)
    assert mean == theta * l0
    assert empirical == predicted
    return empirical, theta, mu4, l0, rb


def build_states(a: np.ndarray, m: int):
    n = len(a)
    q = qvalue(a)
    selectors = list(itertools.combinations(range(n), m))
    states = []
    lookup = {}
    for x in projective_spins(n):
        row = int((a @ x) @ (a @ x))
        raw = int(x @ a @ x)
        for sigma in (-1, 1):
            ls = []
            for s in selectors:
                idx = np.asarray(s)
                child = sigma * int(x[idx] @ a[np.ix_(idx, idx)] @ x[idx])
                ls.append(2 * child - sigma * raw)
            alpha = Fraction(sum(v >= q for v in ls), len(selectors))
            state = dict(x=x.copy(), sigma=sigma, row=row, ls=tuple(ls), alpha=alpha)
            lookup[(sigma, tuple(int(v) for v in x))] = state
            states.append(state)
    return q, selectors, states, lookup


def audit(a: np.ndarray, m: int, name: str):
    n = len(a)
    q, selectors, states, lookup = build_states(a, m)
    active = [d for d in states if d["alpha"]]
    lam = n ** -1.5
    objective = lambda d: -math.log(float(d["alpha"])) + lam * d["row"]
    fmin = min(objective(d) for d in active)
    minimizers = [d for d in active if abs(objective(d) - fmin) < 1e-11]

    pareto = [
        d for d in active
        if not any(
            e["alpha"] >= d["alpha"] and e["row"] <= d["row"]
            and (e["alpha"] > d["alpha"] or e["row"] < d["row"])
            for e in active
        )
    ]
    pareto_types = sorted(set(
        (d["row"], d["alpha"],
         tuple(sorted(v - q for v in d["ls"] if v >= q)))
        for d in pareto
    ), key=str)
    pareto_robust_counts = set()
    pareto_has_empty_every_layer = True
    for d in pareto:
        robust = []
        empty = []
        for k in range(1, n):
            theta, _ = eps_moments(n, k)
            robust.append(sum(theta * value > q for value in d["ls"] if value >= q))
            layer_alphas = []
            for u in itertools.combinations(range(n), k):
                y = d["x"].copy()
                y[list(u)] *= -1
                layer_alphas.append(lookup[(d["sigma"], canonical(y))]["alpha"])
            empty.append(sum(value == 0 for value in layer_alphas))
        pareto_robust_counts.add(tuple(robust))
        pareto_has_empty_every_layer &= all(value > 0 for value in empty)

    profiles = []
    for best in minimizers:
        layer = []
        for k in range(1, n):
            alphas = []
            retained_pairs = 0
            cantelli_sum = Fraction()
            active_indices = [i for i, value in enumerate(best["ls"]) if value >= q]
            for u in itertools.combinations(range(n), k):
                y = best["x"].copy()
                y[list(u)] *= -1
                du = lookup[(best["sigma"], canonical(y))]
                alphas.append(du["alpha"])
                retained_pairs += sum(du["ls"][si] >= q for si in active_indices)

            # Check the exact variance formula and average its Cantelli guarantee
            # over active selectors.  This is arithmetic pair retention only.
            for si in active_indices:
                s = selectors[si]
                idx = np.asarray(s)
                b = -a.copy()
                b[np.ix_(idx, idx)] = a[np.ix_(idx, idx)]
                variance, theta, _mu4, l0, _rb = exact_variance_formula(
                    b, best["x"], best["sigma"], k
                )
                delta = theta * l0 - q
                if delta > 0:
                    cantelli_sum += delta * delta / (variance + delta * delta)
            pair_denom = len(active_indices) * math.comb(n, k)
            actual_pair_retention = Fraction(retained_pairs, pair_denom)
            cantelli_average = cantelli_sum / len(active_indices)
            assert actual_pair_retention >= cantelli_average

            theta, _ = eps_moments(n, k)
            layer.append(
                (k, sum(alpha == 0 for alpha in alphas), len(alphas), str(theta),
                 str(actual_pair_retention), str(cantelli_average))
            )
        profiles.append(tuple(layer))

    # Symmetry-related minimizers have identical layer statistics in these cases.
    distinct_profiles = sorted(set(profiles), key=str)
    slack_hist = {}
    for value in minimizers[0]["ls"]:
        if value >= q:
            slack_hist[value - q] = slack_hist.get(value - q, 0) + 1
    print({
        "name": name,
        "n": n,
        "m": m,
        "q": q,
        "scalar_minimizers": len(minimizers),
        "representative": (
            minimizers[0]["row"],
            int(minimizers[0]["alpha"] * len(selectors)),
            minimizers[0]["sigma"],
            tuple(int(v) for v in minimizers[0]["x"]),
        ),
        "active_slack_hist": slack_hist,
        "pareto_types_(row,alpha,slacks)": pareto_types,
        "pareto_robust_count_profiles": sorted(pareto_robust_counts),
        "all_pareto_columns_have_empty_each_nontrivial_layer": (
            pareto_has_empty_every_layer
        ),
        "distinct_layer_profiles": len(distinct_profiles),
        "representative_layers_(k,empty,total,theta,actual_pair,Cantelli)": profiles[0],
        "empty_count_ranges_by_k": tuple(
            (k, min(profile[k - 1][1] for profile in profiles),
             max(profile[k - 1][1] for profile in profiles), math.comb(n, k))
            for k in range(1, n)
        ),
    })


def main():
    for a, m, name in ((A5, 4, "A5"), (A6, 5, "A6"),
                       (A8, 6, "A8"), (A9, 7, "A9")):
        audit(a, m, name)
    print("PASS noncentral_slack_r46_check")


if __name__ == "__main__":
    main()
