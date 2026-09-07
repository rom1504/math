#!/usr/bin/env python3
"""Exact checks for the Wave 50 switching-circulation obstruction.

The theorem in the accompanying memo is algebraic.  This checker exhausts
the conditional switching mixture, all oriented states, and all relevant
edge marginals on small arbitrary signings and on the stored exact
minimizers.  Every calculation is integral; probabilities are checked after
clearing their power-of-two denominators.
"""

from __future__ import annotations

import itertools
import sys

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from endpoint_transport_r25 import A6, A8, A9  # noqa: E402
from tropical_decomposition_r46_check import (  # noqa: E402
    edges,
    matrix_from_normalized_code,
    qnorm,
    spins,
)


def switch_matrix(A: np.ndarray, U: tuple[int, ...]) -> np.ndarray:
    """Switch A at U, equivalently flip precisely the cut delta(U)."""
    z = np.ones(len(A), dtype=np.int64)
    z[list(U)] = -1
    return (z[:, None] * A * z[None, :]).astype(np.int64)


def cut_vector(n: int, U: tuple[int, ...]) -> np.ndarray:
    """The 0/1 incidence vector of delta(U), in lexicographic edge order."""
    u = set(U)
    return np.asarray([int((i in u) != (j in u)) for i, j in edges(n)], dtype=np.int64)


def oriented_states(A: np.ndarray) -> list[tuple[int, np.ndarray, int]]:
    """Projective x, orientation sigma, and oriented full energy."""
    X = spins(len(A)).astype(np.int64)
    raw = np.einsum("bi,ij,bj->b", X, A.astype(np.int64), X, optimize=True)
    return [(sigma, x, int(sigma * energy)) for sigma in (-1, 1) for x, energy in zip(X, raw)]


def restriction_key(sigma: int, x: np.ndarray, S: tuple[int, ...]) -> tuple[int, ...]:
    """Gauge-free oriented restriction, independent of the signing A."""
    return (sigma,) + tuple(
        int(sigma * x[i] * x[j]) for i, j in itertools.combinations(S, 2)
    )


def conditional_mixture_audit(A: np.ndarray, S: tuple[int, ...]) -> int:
    """Audit the one-selector mixture and Bellman rewrite for every state."""
    n = len(A)
    ee = edges(n)
    q = qnorm(A)
    states = oriented_states(A)
    T = tuple(i for i in range(n) if i not in S)
    subsets = [tuple(T[j] for j in range(len(T)) if (mask >> j) & 1) for mask in range(1 << len(T))]

    # The law is: atom U=empty of mass 1/2, plus a uniform subset of T
    # of total mass 1/2.  Clear the denominator 2^(|T|+1).
    denominator = 1 << (len(T) + 1)
    weights = [1 << len(T)] + [1] * len(subsets)
    outcomes = [tuple()] + subsets
    assert sum(weights) == denominator

    # Conditional edge marginals are 0 on E(S) and 1/4 otherwise.  This
    # separately checks S--T and T--T edges; both have exactly half of the
    # uniform-subset branch toggled.
    marginal_numerators = sum(
        weight * cut_vector(n, U) for weight, U in zip(weights, outcomes)
    )
    expected_numerators = np.asarray(
        [0 if (i in S and j in S) else denominator // 4 for i, j in ee],
        dtype=np.int64,
    )
    assert np.array_equal(marginal_numerators, expected_numerators)

    # Every outcome is a switching symmetry, so its norm is exactly q.
    for U in outcomes:
        assert qnorm(switch_matrix(A, U)) == q

    # Completion slacks s_S are computed exactly by grouping full states by
    # their oriented restriction.
    best_completion_energy: dict[tuple[int, ...], int] = {}
    for sigma, x, energy in states:
        key = restriction_key(sigma, x, S)
        best_completion_energy[key] = max(best_completion_energy.get(key, -10**9), energy)

    selectors = list(itertools.combinations(range(n), len(S)))
    qstar = max(qnorm(A[np.ix_(R, R)]) for R in selectors)

    AS = A[np.ix_(S, S)].astype(np.int64)
    for sigma, x, energy in states:
        slack = q - energy
        child = int(sigma * x[list(S)] @ AS @ x[list(S)])

        # Margin after switching A at U, evaluated at the fixed state d.
        # If F=delta(U), then E_(A^F)(d)=E_A(d^U), so the margin is the
        # original slack at the switched state d^U.
        margin_numerator = 0
        for weight, U in zip(weights, outcomes):
            AU = switch_matrix(A, U)
            switched_energy_at_d = int(sigma * x @ AU @ x)
            z = np.ones(n, dtype=np.int64)
            z[list(U)] = -1
            transformed_energy = int(sigma * (z * x) @ A @ (z * x))
            assert switched_energy_at_d == transformed_energy
            margin_numerator += weight * (q - switched_energy_at_d)

        # E margin = (s(d)+q-c_S(d[S]))/2.  This audits all factors of
        # two from full matrix energy and all factors of four from flips.
        expected_twice = slack + q - child
        assert 2 * margin_numerator == denominator * expected_twice

        key = restriction_key(sigma, x, S)
        sS = q - best_completion_energy[key]
        gS = qstar - child
        jS = sS + gS
        # Bellman rewrite of the same symmetry average.
        bellman_twice = q - qstar + slack - sS + jS
        assert bellman_twice == expected_twice

    # Explicit migrating-ground circulation.  Switching the signing at U and
    # simultaneously switching a fixed parent ground gives a ground again.
    sigma0, x0, energy0 = next(row for row in states if row[2] == q)
    child0 = int(sigma0 * x0[list(S)] @ AS @ x0[list(S)])
    for U in outcomes:
        z = np.ones(n, dtype=np.int64)
        z[list(U)] = -1
        AU = switch_matrix(A, U)
        xU = z * x0
        assert int(sigma0 * xU @ AU @ xU) == q
        # U is disjoint from S, so both the principal signing and the ground
        # restriction on S are unchanged by the circulation.
        assert np.array_equal(AU[np.ix_(S, S)], AS)
        assert np.array_equal(xU[list(S)], x0[list(S)])
        assert int(sigma0 * xU[list(S)] @ AS @ xU[list(S)]) == child0

    return len(states)


def a8_bad_family_audit() -> dict[str, object]:
    """Audit the equal law on the two bad maximal A8 four-selectors."""
    A = A8
    n = len(A)
    ee = edges(n)
    q = qnorm(A)
    states = oriented_states(A)
    bad = [(0, 3, 4, 5), (1, 2, 6, 7)]
    qstar = max(
        qnorm(A[np.ix_(S, S)]) for S in itertools.combinations(range(n), 4)
    )
    assert (q, qstar) == (20, 12)

    # Twice the selector loads w, hence 8p=2-2w.
    two_w = np.zeros(len(ee), dtype=np.int64)
    for S in bad:
        for e in itertools.combinations(S, 2):
            two_w[ee.index(e)] += 1
    eight_p = 2 - two_w
    assert set(map(int, eight_p)) == {1, 2}

    expected_margins = []
    all_support_minima = []
    child_ground_gaps = []
    for sigma, x, energy in states:
        slack = q - energy
        M = np.asarray(
            [int(sigma * A[i, j] * x[i] * x[j]) for i, j in ee],
            dtype=np.int64,
        )
        # Since p=(8p)/8, 4 p.M=(8p).M/2.
        twice_margin = 2 * slack + int(M @ eight_p)
        assert twice_margin % 2 == 0
        expected_margins.append(twice_margin // 2)

    # Every conditional outcome for each bad selector is a switch and has
    # Boolean minimum margin zero.  A ground witness migrates with U and its
    # restriction remains at least four below the maximal child energy.
    for S in bad:
        T = tuple(i for i in range(n) if i not in S)
        AS = A[np.ix_(S, S)].astype(np.int64)
        for mask in range(1 << len(T)):
            U = tuple(T[j] for j in range(len(T)) if (mask >> j) & 1)
            AU = switch_matrix(A, U)
            margins = [q - int(sigma * x @ AU @ x) for sigma, x, _ in states]
            all_support_minima.append(min(margins))
            assert min(margins) == 0

            sigma0, x0, _ = next(row for row in states if row[2] == q)
            z = np.ones(n, dtype=np.int64)
            z[list(U)] = -1
            xU = z * x0
            assert int(sigma0 * xU @ AU @ xU) == q
            child = int(sigma0 * xU[list(S)] @ AS @ xU[list(S)])
            child_ground_gaps.append(qstar - child)
            assert qstar - child >= 4

    assert min(expected_margins) == 8
    return {
        "q": q,
        "qstar": qstar,
        "bad_selectors": bad,
        "eight_p_values": sorted(set(map(int, eight_p))),
        "fractional_minimum_margin": min(expected_margins),
        "switching_support_outcomes_checked": len(all_support_minima),
        "support_boolean_minima": sorted(set(all_support_minima)),
        "migrated_ground_child_gap_values": sorted(set(child_ground_gaps)),
    }


def main() -> None:
    checks = 0

    # Exhaust every switching-normalized signing through order four and every
    # nonempty proper selector.  This is a generic identity audit, not a
    # minimizer-only computation.
    arbitrary_signings = 0
    for n in range(2, 5):
        codes = 1 << ((n - 1) * (n - 2) // 2)
        for code in range(codes):
            A = matrix_from_normalized_code(n, code)
            arbitrary_signings += 1
            for m in range(1, n):
                for S in itertools.combinations(range(n), m):
                    checks += conditional_mixture_audit(A, S)

    # Larger stored exact minimizers: audit representative small, middle, and
    # near-diagonal selectors so both one-outside and two-outside edge cases
    # occur at nontrivial sizes.
    stored_cases = [
        (A6, (0, 1, 2)),
        (A8, (0, 3, 4, 5)),
        (A8, (0, 1, 2, 5, 6, 7)),
        (A9, (0, 1, 2, 3, 4, 5, 6, 7)),
    ]
    for A, S in stored_cases:
        checks += conditional_mixture_audit(A, S)

    a8 = a8_bad_family_audit()
    print("arbitrary normalized signings audited", arbitrary_signings)
    print("oriented state/selector identities checked", checks)
    print("A8 bad-family switching circulation", a8)
    print("PASS global_all_bad_r50_check")


if __name__ == "__main__":
    main()
