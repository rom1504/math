#!/usr/bin/env python3
"""Wave 42 audit of the fixed-density reveal/affinity identity."""

from __future__ import annotations

import math
import sys
from collections import defaultdict
from fractions import Fraction

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from endpoint_transport_r25 import A6, A8, A9, restriction_key
from harmonic_resonance_r34_check import Endpoint


def audit(matrix: np.ndarray, beta: float, m: int) -> dict[str, float | int]:
    Hm = Endpoint(matrix, beta, m)
    Hp = Endpoint(matrix, beta, m + 1)
    bm = np.exp(-Hm.Fall)
    bp = np.exp(-Hp.Fall)
    map_m = {S: j for j, S in enumerate(Hm.selectors)}
    map_p = {R: j for j, R in enumerate(Hp.selectors)}
    Zm = np.sum(bm, axis=0)

    identity_error = 0.0
    composition_error = 0.0
    crossing_slack_error = 0.0
    crossing_count = 0
    maximum_crossing_odds = 0.0
    maximum_crossing_jump = 0.0
    maximum_reveal_cost = 0.0
    maximum_affinity_cost = 0.0
    maximum_cross_level_cost = 0.0
    minimum_affinity = 1.0
    records: list[tuple[float, ...]] = []

    for i in range(1, Hm.n):
        containing = [R for R in Hp.selectors if i in R]
        for x, y in Hm.edges[i]:
            numerator = 0.0
            geometric_terms = []
            reveal_factors = []
            zpx = 0.0
            zpy = 0.0
            for R in containing:
                S = tuple(j for j in R if j != i)
                js = map_m[S]
                jr = map_p[R]
                b_sx = float(bm[js, x])
                b_sy = float(bm[js, y])
                b_rx = float(bp[jr, x])
                b_ry = float(bp[jr, y])
                composition_error = max(composition_error, abs(b_sx - b_sy))
                geom = math.sqrt(b_rx * b_ry)
                reveal = b_sx / geom
                assert reveal <= 1.0 + 2e-11
                numerator += b_sx
                geometric_terms.append(geom)
                reveal_factors.append(reveal)
                zpx += b_rx
                zpy += b_ry

            geom_a = np.asarray(geometric_terms)
            reveal_a = np.asarray(reveal_factors)
            affinity = float(np.sum(geom_a) / math.sqrt(zpx * zpy))
            reveal_mean = float(np.sum(geom_a * reveal_a) / np.sum(geom_a))
            cross_level = math.sqrt(zpx * zpy / (float(Zm[x]) * float(Zm[y])))
            rx = numerator / float(Zm[x])
            ry = numerator / float(Zm[y])
            lhs = math.sqrt(rx * ry)
            rhs = cross_level * affinity * reveal_mean
            identity_error = max(identity_error, abs(lhs - rhs))

            reveal_cost = -math.log(reveal_mean)
            affinity_cost = -math.log(affinity)
            cross_level_cost = -math.log(cross_level)
            total_cost = reveal_cost + affinity_cost + cross_level_cost
            assert abs(total_cost + math.log(lhs)) < 3e-10
            minimum_affinity = min(minimum_affinity, affinity)
            maximum_reveal_cost = max(maximum_reveal_cost, reveal_cost)
            maximum_affinity_cost = max(maximum_affinity_cost, affinity_cost)
            maximum_cross_level_cost = max(maximum_cross_level_cost, cross_level_cost)

            omega = float(Hm.lognu[y] - Hm.lognu[x])
            chi = float(Hm.g[y] - Hm.g[x])
            crossing = omega * (omega + chi) <= 1e-13 and abs(chi) > 1e-12
            if crossing:
                crossing_count += 1
                maximum_crossing_odds = max(maximum_crossing_odds, abs(omega))
                maximum_crossing_jump = max(maximum_crossing_jump, abs(chi))
                crossing_slack_error = max(
                    crossing_slack_error, abs(omega) / 2.0 - total_cost
                )
                records.append(
                    (
                        abs(omega),
                        abs(chi),
                        total_cost,
                        reveal_cost,
                        affinity_cost,
                        cross_level_cost,
                    )
                )

    assert identity_error < 3e-11
    assert composition_error < 3e-11
    assert crossing_slack_error < 3e-10
    records.sort(reverse=True)
    largest = records[0] if records else (0.0,) * 6
    return {
        "n": Hm.n,
        "m": m,
        "beta": beta,
        "identity_error": identity_error,
        "composition_invariance_error": composition_error,
        "crossing_slack_error": crossing_slack_error,
        "crossing_count": crossing_count,
        "maximum_crossing_odds": maximum_crossing_odds,
        "maximum_crossing_jump": maximum_crossing_jump,
        "minimum_affinity": minimum_affinity,
        "maximum_reveal_cost": maximum_reveal_cost,
        "maximum_affinity_cost": maximum_affinity_cost,
        "maximum_cross_level_cost": maximum_cross_level_cost,
        "largest_crossing_record": largest,
    }


def tropical_component(H: Endpoint, state: int) -> tuple[int, Fraction]:
    """Leading low-temperature exponent and coefficient of U at one state."""
    rows: list[tuple[int, int]] = []
    for S in H.selectors:
        key = restriction_key(H.cuts[state], S)
        inds = [
            j
            for j, d in enumerate(H.cuts)
            if restriction_key(d, S) == key
        ]
        child = int(
            sum(
                H.A[a, b] * H.cuts[state, a, b]
                for a in S
                for b in S
            )
        )
        external = [int(H.energies[j] - child) for j in inds]
        maximum = max(external)
        rows.append((maximum, external.count(maximum)))
    exponent = min(maximum for maximum, _ in rows)
    coefficient = sum(
        (Fraction(1, multiplicity) for maximum, multiplicity in rows if maximum == exponent),
        Fraction(0, 1),
    )
    return exponent, coefficient


def tropical_exponents(H: Endpoint) -> np.ndarray:
    """Compute min_S max external energy simultaneously for every state."""
    exponent = np.full(H.N, 10**9, dtype=np.int64)
    for S in H.selectors:
        groups: dict[tuple[int, ...], list[int]] = defaultdict(list)
        for j, d in enumerate(H.cuts):
            groups[restriction_key(d, S)].append(j)
        for inds in groups.values():
            first = inds[0]
            child = int(
                sum(
                    H.A[a, b] * H.cuts[first, a, b]
                    for a in S
                    for b in S
                )
            )
            maximum = max(int(H.energies[j] - child) for j in inds)
            exponent[inds] = np.minimum(exponent[inds], maximum)
    return exponent


def exact_fixed_density_wall() -> dict[str, object]:
    """A9,m=4 zero-base-odds edge with an unbounded low-temperature jump."""
    H = Endpoint(A9, 2.0, 4)
    x, y, bit = 292, 308, 4
    assert int(np.max(np.abs(H.energies))) == 24
    assert int(H.energies[x]) == int(H.energies[y]) == -8
    assert (x, y) in map(tuple, H.edges[bit])
    ax, cx = tropical_component(H, x)
    ay, cy = tropical_component(H, y)
    assert (ax, ay) == (16, 8)
    assert cx == Fraction(1318, 105)
    assert cy == Fraction(1, 8)
    ratio = cy / cx
    assert ratio == Fraction(105, 10544)
    beta = 8.0
    Hb = Endpoint(A9, beta, 4)
    jump = float(Hb.g[y] - Hb.g[x])
    asymptotic = 8.0 * beta + math.log(float(ratio))
    assert abs(jump - asymptotic) < 3e-9
    exponents = tropical_exponents(H)
    grounds = np.flatnonzero(H.energies == 24)
    assert set(map(int, exponents[grounds])) == {12}
    assert int(np.max(H.energies - exponents)) == 12
    # The ground line 24-12s dominates for every 0<=s<=1.  The y endpoint
    # has line -8-8s, whose smallest gap is 28 at s=1; x is still rarer.
    rarity_exponent = min(
        (24 - 12 * s) - (-8 - ay * s) for s in (0, 1)
    )
    assert rarity_exponent == 28
    return {
        "states": (x, y),
        "bit": bit,
        "parent_energies": (int(H.energies[x]), int(H.energies[y])),
        "base_odds_all_beta": 0,
        "tropical_exponents": (ax, ay),
        "tropical_coefficients": (str(cx), str(cy)),
        "coefficient_ratio": str(ratio),
        "score_jump_asymptotic": "8 beta + log(105/10544) + o(1)",
        "beta8_jump": jump,
        "beta8_asymptotic": asymptotic,
        "uniform_context_rarity_exponent": rarity_exponent,
        "edge_energy_bound": "O(beta^2 exp(-28 beta))",
    }


def main() -> None:
    cases = [
        ("A6", A6, 0.5, 3),
        ("A6", A6, 2.0, 3),
        ("A6", A6, 2.0, 4),
        ("A8", A8, 0.5, 4),
        ("A8", A8, 2.0, 4),
        ("A8", A8, 2.0, 6),
        ("A9", A9, 0.5, 4),
        ("A9", A9, 2.0, 4),
        ("A9", A9, 2.0, 7),
    ]
    for name, matrix, beta, m in cases:
        print(name, audit(matrix, beta, m))
    print("exact fixed-density wall", exact_fixed_density_wall())
    print("PASS harmonic_fixed_density_reveal_r42_check")


if __name__ == "__main__":
    main()
