#!/usr/bin/env python3
"""Exact checks for the Wave 46 common-core escape calculation."""

from __future__ import annotations

import math
import sys
from fractions import Fraction

import numpy as np

sys.path.insert(0, "/home/math/quadra/tmp")
from envelope_block_cover_r27 import A6, A8, A9
from hard_center_branching_r45_check import audit, down_up_kernel
from high_replica_pressure_r44_explore import geometry


def finite_audit(name: str, matrix: np.ndarray, m: int) -> None:
    old = audit(name, matrix, m)
    j2_over_j1 = old["second_ratio"]
    max_degree = old["max_degree"]
    print(f"{name}: J2/J1={j2_over_j1}, max-a={max_degree}")
    normalized = []
    for ell, rho, lam, _ in old["rows"]:
        delta = 1 - lam
        escape = 1 - rho
        basic = (rho - lam) / delta
        # From the complete Johnson spectrum:
        # Xi/J1 = delta J2/J1 - (rho-lambda).
        xi_over_j1 = delta * j2_over_j1 - (rho - lam)
        assert xi_over_j1 >= 0
        assert basic + xi_over_j1 / delta == j2_over_j1
        assert max_degree >= j2_over_j1 >= basic
        normalized.append((ell, rho, lam, escape, basic, xi_over_j1))
        print(
            f"  ell={ell}: escape={escape}, delta={delta}, "
            f"basic={basic}, Xi/J1={xi_over_j1}"
        )

    # A nonnegative kernel mixture cannot beat its best constituent's
    # normalized excess.  Use deterministic positive rational weights.
    raw_weights = [Fraction((i + 1) ** 2, 1) for i in range(len(normalized))]
    total = sum(raw_weights)
    weights = [w / total for w in raw_weights]
    rho_mix = sum(w * row[1] for w, row in zip(weights, normalized))
    lam_mix = sum(w * row[2] for w, row in zip(weights, normalized))
    delta_mix = 1 - lam_mix
    basic_mix = (rho_mix - lam_mix) / delta_mix
    weighted_basic = sum(
        w * (1 - row[2]) * row[4] for w, row in zip(weights, normalized)
    ) / delta_mix
    assert basic_mix == weighted_basic
    assert basic_mix <= max(row[4] for row in normalized)

    # Exact layer-cake/coarea identity at cap H=4.  These examples have
    # deficits in 4 Z, so the positive-slack levels 0<h<4 all select D=0.
    z, selector_tuples, _, _, row_cost, deficits, _ = geometry(matrix, m)
    n = len(matrix)
    eligible = row_cost <= 2 * n * (n - 1)
    assert np.all(eligible)
    selectors = list(map(frozenset, selector_tuples))
    cap = 4
    clipped = np.minimum(deficits.astype(int), cap)
    slack = cap - clipped
    mean_slack = Fraction(int(np.sum(slack)), slack.size)
    for ell, rho, lam, escape, _, _ in normalized:
        kernel = down_up_kernel(selectors, n, m, ell)
        abs_sum = Fraction(0)
        directed_layer_sum = Fraction(0)
        for iz in range(len(z)):
            for i in range(len(selectors)):
                for j in range(len(selectors)):
                    kij = kernel[i, j]
                    abs_sum += kij * abs(int(clipped[iz, i]) - int(clipped[iz, j]))
                    # Integral from 0 to cap of 1{D_S<=h<D_T} dh.
                    directed_layer_sum += kij * max(
                        0, int(clipped[iz, j]) - int(clipped[iz, i])
                    )
        norm = len(z) * len(selectors)
        half_l1 = abs_sum / (2 * norm)
        directed = directed_layer_sum / norm
        assert half_l1 == directed
        # Since the cap-four layer is exactly the D=0 hard family, the
        # integrated escape divided by integrated incidence is 1-rho.
        assert directed / mean_slack == escape


def falling(x: int, r: int) -> int:
    out = 1
    for j in range(r):
        out *= x - j
    return out


def factorial_dirichlet_comparison() -> None:
    """Check 1-lambda_j(m-s) <= s(1-lambda_j(m-1)) exactly."""
    for n in range(4, 31):
        for m in range(2, n):
            codim = n - m
            for s in range(1, m + 1):
                for j in range(1, min(m, codim) + 1):
                    lam_s = Fraction(
                        falling(m - s, j) * falling(codim, j),
                        falling(m, j) * falling(codim + s, j),
                    )
                    lam_1 = Fraction(
                        falling(m - 1, j) * falling(codim, j),
                        falling(m, j) * falling(codim + 1, j),
                    )
                    assert 1 - lam_s <= s * (1 - lam_1)


def cylinder_retention(n: int, m: int, s: int, r: int) -> Fraction:
    """Retention of F={S: W subset S}, |W|=r, under K_{m-s}."""
    total = Fraction(0)
    available = n - m + s
    for j in range(min(r, s) + 1):
        if s - j > m - r:
            continue
        removed = Fraction(
            math.comb(r, j) * math.comb(m - r, s - j), math.comb(m, s)
        )
        reselect = Fraction(falling(s, j), falling(available, j))
        total += removed * reselect
    return total


def scalable_cylinder_audit() -> None:
    print("scalable fixed-coordinate cylinder")
    previous_retention = None
    for scale in (8, 12, 16, 20):
        n = scale**4
        m = n // 2
        s = n // scale
        r = scale
        density = Fraction(falling(m, r), falling(n, r))
        retention = cylinder_retention(n, m, s, r)
        delta = Fraction(s * n, m * (n - m + s))
        lam = 1 - delta
        basic = (retention - lam) / delta
        assert retention >= density  # positive correlation from the common core
        assert basic <= density       # full-spectrum correction is indispensable
        if scale >= 12:
            assert retention < lam    # near-core spectral excess is negative
        if previous_retention is not None:
            # Numerical monotonicity is not used in the proof, but catches
            # mistakes in the exact hypergeometric formula.
            assert float(retention) < float(previous_retention)
        previous_retention = retention
        print(
            f"  L={scale}: -log(a)/L={-math.log(float(density))/scale:.8f}, "
            f"retention={float(retention):.8f}, escape={float(1-retention):.8f}, "
            f"L*delta={float(scale*delta):.8f}, normalized={float(basic):.8f}"
        )


def main() -> None:
    factorial_dirichlet_comparison()
    for name, matrix, m in (("A6", A6, 5), ("A8", A8, 6), ("A9", A9, 7)):
        finite_audit(name, matrix, m)
    scalable_cylinder_audit()
    print("PASS: escape identity, normalized/full-spectrum extraction, mixture collapse, coarea, and cylinder model")


if __name__ == "__main__":
    main()
