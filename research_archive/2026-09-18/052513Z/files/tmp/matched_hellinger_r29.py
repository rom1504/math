#!/usr/bin/env python3
"""Exact matched-temperature endpoint limits on finite signing examples.

All temporary output belongs to Wave 29 Route 3.  The calculation uses
oriented cuts, matching the normalizations in endpoint_transport_r25.py.
"""

from __future__ import annotations

from fractions import Fraction
import math
import sys

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from endpoint_transport_r25 import A4, A6, A8, A9, oriented_cuts, restriction_key
from check_finite_bridge_r16 import A5
from parent_cut_entropy_r26 import Endpoint


def matched_zero_temperature(a_mat: np.ndarray) -> dict[str, object]:
    """Return the exact beta=gamma -> infinity endpoint law.

    For one deletion, first retain only selectors with maximal child optimum;
    then choose uniformly among all child-ground pairs on those selectors and
    uniformly among the best parent extensions of the chosen child cut.
    """
    n = len(a_mat)
    cuts = oriented_cuts(n)
    nd = len(cuts)
    energy = np.einsum("bij,ij->b", cuts, a_mat, optimize=True).astype(int)
    selectors = [tuple(j for j in range(n) if j != i) for i in range(n)]

    selector_data = []
    for ids in selectors:
        child_energy = np.asarray([
            int(np.sum(a_mat[np.ix_(ids, ids)] * d[np.ix_(ids, ids)]))
            for d in cuts
        ])
        groups: dict[tuple[int, ...], list[int]] = {}
        for di, d in enumerate(cuts):
            groups.setdefault(restriction_key(d, ids), []).append(di)
        group_list = list(groups.values())
        child_values = [int(child_energy[g[0]]) for g in group_list]
        qs = max(child_values)
        child_grounds = [g for g, value in zip(group_list, child_values) if value == qs]
        selector_data.append((qs, child_grounds))

    q_star = max(qs for qs, _ in selector_data)
    active = [(si, grounds) for si, (qs, grounds) in enumerate(selector_data) if qs == q_star]
    pair_count = sum(len(grounds) for _, grounds in active)

    output = [Fraction(0) for _ in range(nd)]
    channel = [[Fraction(0) for _ in range(nd)] for _ in range(n)]
    weights = [Fraction(0) for _ in range(n)]
    extension_profile: dict[tuple[int, int], int] = {}
    for si, grounds in active:
        weights[si] = Fraction(n * len(grounds), pair_count)
        for group in grounds:
            best = max(int(energy[di]) for di in group)
            maximizers = [di for di in group if int(energy[di]) == best]
            extension_profile[(best, len(maximizers))] = extension_profile.get((best, len(maximizers)), 0) + 1
            for di in maximizers:
                mass = Fraction(1, pair_count * len(maximizers))
                output[di] += mass
                # q_S itself chooses 1/g_S, then 1/k extensions.
                channel[si][di] += Fraction(1, len(grounds) * len(maximizers))

    assert sum(output) == 1
    for si, _ in active:
        assert sum(channel[si]) == 1
    assert sum(weights, Fraction(0)) == n

    q_parent = max(int(e) for e in energy)
    ground = [di for di, e in enumerate(energy) if int(e) == q_parent]
    ground_count = len(ground)
    ebar = sum((output[di] * int(energy[di]) for di in range(nd)), Fraction(0))
    slope = Fraction(q_parent) - ebar

    entropy = -sum(float(p) * math.log(float(p)) for p in output if p)
    constant = math.log(ground_count) - entropy
    if slope == 0:
        assert all(output[di] == 0 for di in range(nd) if di not in ground)
        assert constant >= -2e-14

    directed = [(s, t) for s in range(n) for t in range(n) if s != t]
    finite_h = 0.0
    for s, t in directed:
        mu_s = [weights[s] * p for p in channel[s]]
        mu_t = [weights[t] * p for p in channel[t]]
        finite_h += 0.5 * sum(
            (math.sqrt(float(x)) - math.sqrt(float(y))) ** 2
            for x, y in zip(mu_s, mu_t)
        ) / len(directed)

    positive_masses: dict[Fraction, int] = {}
    for p in output:
        if p:
            positive_masses[p] = positive_masses.get(p, 0) + 1

    return {
        "q_parent": q_parent,
        "q_star": q_star,
        "active_selectors": len(active),
        "child_ground_pairs": pair_count,
        "weights": tuple(weights),
        "extension_profile": extension_profile,
        "ground_count": ground_count,
        "output_support": sum(bool(p) for p in output),
        "output_mass_profile": positive_masses,
        "ebar": ebar,
        "slope": slope,
        "constant": constant,
        "mean_finite_h": finite_h,
    }


def spin_from_oriented_cut(d: np.ndarray) -> np.ndarray:
    """Recover x_0=1 from d=sigma xx^T (off diagonal), n>=3."""
    sigma = int(d[0, 1] * d[1, 2] * d[2, 0])
    assert sigma in (-1, 1)
    x = np.ones(len(d), dtype=np.int64)
    x[1:] = d[0, 1:] // sigma
    check = sigma * np.outer(x, x)
    np.fill_diagonal(check, 0)
    assert np.array_equal(check, d)
    return x


def matched_high_temperature(a_mat: np.ndarray) -> dict[str, object]:
    """Exact beta^4 coefficients for beta=gamma near zero."""
    n = len(a_mat)
    m = n - 1
    k = n - m
    cuts = oriented_cuts(n)
    selectors = [tuple(j for j in range(n) if j != i) for i in range(n)]
    gs = np.empty((len(selectors), len(cuts)), dtype=np.int64)
    for di, d in enumerate(cuts):
        x = spin_from_oriented_cut(d)
        for si, ids in enumerate(selectors):
            outside = tuple(j for j in range(n) if j not in ids)
            boundary_square = sum(int(np.dot(a_mat[u, list(ids)], x[list(ids)])) ** 2 for u in outside)
            gs[si, di] = m * k - boundary_square

    gbar = np.mean(gs, axis=0)
    row_square = np.asarray([
        int(np.sum((a_mat * d).sum(axis=1) ** 2)) for d in cuts
    ])
    theta = Fraction(k * m * (m - 1), n * (n - 1) * (n - 2))
    predicted_gbar = float(theta) * (n * (n - 1) - row_square)
    assert np.max(np.abs(gbar - predicted_gbar)) < 1e-12
    assert abs(float(np.mean(gbar))) < 1e-12

    parent_beta4 = 2.0 * float(np.mean(gbar**2))
    directed = [(s, t) for s in range(n) for t in range(n) if s != t]
    hell_beta4 = 0.5 * float(np.mean([
        np.mean((gs[s] - gs[t]) ** 2) for s, t in directed
    ]))

    # Check the second derivative of the exact endpoint likelihood directly.
    beta = 1e-4
    endpoint = Endpoint(a_mat).evaluate(beta, beta)
    second = (np.asarray(endpoint["f"]) - 1.0) / beta**2
    assert np.max(np.abs(second - 2.0 * gbar)) < 3e-2

    return {
        "theta": theta,
        "row_square_profile": {int(v): int(np.sum(row_square == v)) for v in np.unique(row_square)},
        "parent_entropy_beta4_coefficient": parent_beta4,
        "selector_hellinger_beta4_coefficient": hell_beta4,
    }


def main() -> None:
    for name, a_mat in (("A4", A4), ("A5", A5), ("A6", A6), ("A8", A8), ("A9", A9)):
        z = matched_zero_temperature(a_mat)
        print(name)
        for key, value in z.items():
            print(" ", key, value)
        assert z["slope"] == 0
        # At beta=8 the exact finite endpoint is already close to the limit.
        beta8 = Endpoint(a_mat).evaluate(8.0, 8.0)
        assert abs(float(beta8["ent_cut"]) - float(z["constant"])) < 2e-5
        assert abs(float(beta8["selector_h_mean"]) - float(z["mean_finite_h"])) < 5e-4
        high = matched_high_temperature(a_mat)
        print(" high-temperature")
        for key, value in high.items():
            print("  ", key, value)
    print("PASS matched zero-temperature formula and finite audits")


if __name__ == "__main__":
    main()
