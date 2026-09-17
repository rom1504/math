#!/usr/bin/env python3
"""Exact checks for the Wave 49 full-slack cross-block calculation.

All arithmetic in the profile and margin checks is integral (or doubled
integer arithmetic for the 1/2 selector weights).  The only LP fact used in
the accompanying memo also has explicit primal and dual certificates here.
"""

from __future__ import annotations

import itertools
import sys
from collections import Counter

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from endpoint_transport_r25 import A8  # noqa: E402
from tropical_decomposition_r46_check import edges, qnorm, spins  # noqa: E402


def full_states(A: np.ndarray) -> tuple[list[tuple[int, ...]], np.ndarray, np.ndarray, int]:
    """Return the distinct oriented edge-feature states, energies, and slacks."""
    n = len(A)
    ee = edges(n)
    X = spins(n).astype(np.int64)
    raw = np.einsum("bi,ij,bj->b", X, A.astype(np.int64), X, optimize=True)
    q = qnorm(A)
    state_by_feature: dict[tuple[int, ...], int] = {}
    for sigma in (-1, 1):
        for x, raw_energy in zip(X, raw):
            feature = tuple(
                int(sigma * A[i, j] * x[i] * x[j]) for i, j in ee
            )
            energy = int(sigma * raw_energy)
            if feature in state_by_feature:
                assert state_by_feature[feature] == energy
            state_by_feature[feature] = energy
    features = list(state_by_feature)
    M = np.asarray(features, dtype=np.int64)
    energy = np.asarray([state_by_feature[f] for f in features], dtype=np.int64)
    slack = q - energy
    assert np.all(2 * np.sum(M, axis=1) == energy)
    return features, M, slack, q


def restriction_key(feature: tuple[int, ...], ids: list[int]) -> tuple[int, ...]:
    return tuple(feature[i] for i in ids)


def selector_profile(
    A: np.ndarray,
    features: list[tuple[int, ...]],
    slack: np.ndarray,
    S: tuple[int, ...],
    qstar: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Evaluate s_S, g_S, and j_S=s_S+g_S on every full state."""
    ee = edges(len(A))
    internal = list(itertools.combinations(S, 2))
    ids = [ee.index(e) for e in internal]

    # Bellman completion slack: minimum parent slack over each restriction.
    best: dict[tuple[int, ...], int] = {}
    for feature, s in zip(features, slack):
        key = restriction_key(feature, ids)
        best[key] = min(best.get(key, 10**9), int(s))

    sS = np.asarray([best[restriction_key(f, ids)] for f in features], dtype=np.int64)
    cS = np.asarray(
        [2 * sum(f[i] for i in ids) for f in features], dtype=np.int64
    )
    gS = qstar - cS
    jS = sS + gS
    assert np.all(sS >= 0)
    assert np.all(gS >= 0)
    assert np.all(sS <= slack)  # Bellman projection inequality.
    assert np.all(jS % 4 == 0)
    return sS, gS, jS


def weighted_bellman_identity() -> dict[str, object]:
    """Check the exact selector-load identity on the two bad A8 blocks."""
    A = A8
    ee = edges(len(A))
    features, M, slack, q = full_states(A)
    selectors = [(0, 3, 4, 5), (1, 2, 6, 7)]
    qstar = max(
        qnorm(A[np.ix_(S, S)])
        for S in itertools.combinations(range(len(A)), 4)
    )
    assert (q, qstar) == (20, 12)

    profiles = [selector_profile(A, features, slack, S, qstar) for S in selectors]
    deltas = [int(np.min(profile[2])) for profile in profiles]
    assert deltas == [4, 4]

    # lambda=(1/2,1/2).  Store 2w integrally.
    two_w = np.zeros(len(ee), dtype=np.int64)
    for S in selectors:
        for e in itertools.combinations(S, 2):
            two_w[ee.index(e)] += 1
    # p=(1-w)/4, so 8p=2-2w is integral.
    eight_p = 2 - two_w
    assert set(map(int, eight_p)) == {1, 2}

    # Twice the lambda-averages are integral.
    two_sbar = profiles[0][0] + profiles[1][0]
    two_jbar = profiles[0][2] + profiles[1][2]

    # 2G = 2s + sum_e (8p_e) M_e, because G=s+4 p.M.
    two_G = 2 * slack + M @ eight_p
    # Twice the exact identity
    # G=(q-qstar+jbar+s-sbar)/2 is
    # 2G=q-qstar+jbar+s-sbar.  Multiply once more by 2.
    rhs_times_two = 2 * (q - qstar + slack) + two_jbar - two_sbar
    assert np.array_equal(2 * two_G, rhs_times_two)
    assert np.all(two_G % 2 == 0)
    G = two_G // 2
    assert int(np.min(G)) == 8

    union = {
        e for S in selectors for e in itertools.combinations(S, 2)
    }
    threshold_sets = {
        "empty": [],
        "outside_union": [e for e in ee if e not in union],
        "all_edges": ee,
    }
    threshold_results: dict[str, dict[str, object]] = {}
    for name, F in threshold_sets.items():
        ids = [ee.index(e) for e in F]
        margin = slack + (4 * np.sum(M[:, ids], axis=1) if ids else 0)
        minimum = int(np.min(margin))
        witness_slacks = Counter(int(s) for s, value in zip(slack, margin) if value == minimum)
        assert minimum == 0
        threshold_results[name] = {
            "edge_count": len(F),
            "minimum_integral_margin": minimum,
            "witness_slacks": dict(sorted(witness_slacks.items())),
        }

    return {
        "q": q,
        "qstar": qstar,
        "bad_selectors": selectors,
        "joint_gap_deltas": deltas,
        "canonical_p_values": sorted({f"{v}/8" for v in map(int, eight_p)}),
        "canonical_fractional_minimum_margin": int(np.min(G)),
        "canonical_minimizer_count": int(np.sum(G == np.min(G))),
        "theorem_lower_bound": (q - qstar + min(deltas)) // 2,
        "threshold_roundings": threshold_results,
    }


def pressure_game() -> dict[str, object]:
    """Compute the exact full-slack fractional/Boolean gap on four A8 edges."""
    A = A8
    ee = edges(len(A))
    features, M, slack, q = full_states(A)
    H = [(0, 5), (3, 4), (1, 6), (2, 7)]
    ids = [ee.index(e) for e in H]
    MH = M[:, ids]

    # Primal certificate p_e=1/2: every state has margin at least four.
    half_margin = slack + 2 * np.sum(MH, axis=1)
    assert int(np.min(half_margin)) == 4

    # Dual certificate: a half/half law on an actual ground and its displayed
    # slack-eight H-antipode has E[s]=4 and zero edge means.  Hence the exact
    # fractional value is at most four, matching the primal certificate.
    dual_pair: tuple[int, int] | None = None
    for i in np.flatnonzero(slack == 0):
        for j in np.flatnonzero(slack == 8):
            if np.array_equal(MH[j], -MH[i]):
                dual_pair = (int(i), int(j))
                break
        if dual_pair is not None:
            break
    assert dual_pair is not None
    i, j = dual_pair
    assert np.array_equal(MH[i] + MH[j], np.zeros(len(H), dtype=np.int64))
    assert int(slack[i] + slack[j]) == 8

    # Exhaust the 16 Boolean choices.  Exact minimality predicts optimum zero;
    # this small restricted game checks it directly without using that fact.
    boolean_values: dict[int, int] = {}
    best = -10**9
    best_masks: list[int] = []
    for mask in range(1 << len(H)):
        z = np.asarray([(mask >> k) & 1 for k in range(len(H))], dtype=np.int64)
        value = int(np.min(slack + 4 * (MH @ z)))
        boolean_values[value] = boolean_values.get(value, 0) + 1
        if value > best:
            best = value
            best_masks = [mask]
        elif value == best:
            best_masks.append(mask)
    assert best == 0

    layer_table = Counter(
        (int(s), int(np.sum(row)), int(value))
        for s, row, value in zip(slack, MH, half_margin)
    )
    return {
        "pressure_edges": H,
        "fractional_primal": "p_e=1/2",
        "fractional_minimum_margin": int(np.min(half_margin)),
        "dual_slacks": [int(slack[i]), int(slack[j])],
        "dual_H_features": [list(map(int, MH[i])), list(map(int, MH[j]))],
        "boolean_maximum_margin": best,
        "boolean_best_masks": best_masks,
        "boolean_margin_histogram": dict(sorted(boolean_values.items())),
        "half_flip_layer_table": {
            str(key): count for key, count in sorted(layer_table.items())
        },
    }


def main() -> None:
    print("WEIGHTED BELLMAN FULL-SLACK IDENTITY")
    print(weighted_bellman_identity())
    print("FOUR-EDGE FULL-SLACK INTEGRALITY GAP")
    print(pressure_game())
    print("PASS")


if __name__ == "__main__":
    main()
