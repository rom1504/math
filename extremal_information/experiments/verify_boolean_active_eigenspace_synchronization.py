#!/usr/bin/env python3
"""Finite checks for boolean_active_eigenspace_synchronization.md."""

from __future__ import annotations

from itertools import product
from math import acos, cos, pi, sqrt

import numpy as np


def verify_circle_covering() -> None:
    for count in (8, 16, 32, 64):
        # Equally spaced oriented caps have covering angle pi/count.
        epsilon = 1 - cos(pi / count)
        lower = (2 * epsilon) ** -0.5
        assert count + 1e-12 >= lower
        assert count * sqrt(epsilon) > 2
    print("sphere-covering exponent on S^1: PASS")


def verify_trust_gap_decomposition() -> None:
    n = 4
    # V consists of vectors constant on each coordinate pair.
    v1 = np.asarray((1, 1, 0, 0), dtype=float) / sqrt(2)
    v2 = np.asarray((0, 0, 1, 1), dtype=float) / sqrt(2)
    basis = np.stack((v1, v2), axis=1)
    projection = basis @ basis.T
    spins = np.asarray(list(product((-1, 1), repeat=n)), dtype=float)
    q_values = np.einsum("bi,ij,bj->b", spins, projection, spins) / n
    projected = spins @ basis / sqrt(n)
    beta = 0.7

    worst_gap = 0.0
    worst_support = 1.0
    for theta in np.linspace(0, 2 * pi, 4097)[:-1]:
        direction = np.asarray((cos(theta), np.sin(theta)))
        supports = projected @ direction
        values = q_values / 2 + beta * supports
        index = int(np.argmax(values))
        gap = 0.5 + beta - float(values[index])
        defect = 1 - float(q_values[index])
        support_deficit = 1 - float(supports[index])
        assert abs(gap - (defect / 2 + beta * support_deficit)) < 1e-12
        assert defect <= 2 * gap + 1e-12
        assert support_deficit <= gap / beta + 1e-12
        worst_gap = max(worst_gap, gap)
        worst_support = min(worst_support, float(supports[index]))

    exact_witnesses = spins[np.abs(q_values - 1) < 1e-12]
    assert len(exact_witnesses) == 4 == 2**2
    epsilon = 1 - 1 / sqrt(2)
    assert worst_gap <= beta * epsilon + 1e-9
    assert abs(worst_support - 1 / sqrt(2)) < 2e-3
    print(
        "active trust decomposition:",
        f"Gamma={worst_gap:.6f}, exact_witnesses={len(exact_witnesses)}",
    )


def verify_finite_extraction_geometry() -> None:
    epsilon = 1e-3
    rho = sqrt(epsilon)
    # A source witness within support deficit epsilon of a net centre.
    centre_angle = 0.37
    witness_angle = centre_angle + acos(1 - epsilon)
    u = np.asarray((cos(witness_angle), np.sin(witness_angle)))
    centre = np.asarray((cos(centre_angle), np.sin(centre_angle)))
    assert centre @ u >= 1 - epsilon - 1e-12

    target = np.asarray((cos(centre_angle + rho), np.sin(centre_angle + rho)))
    deficit = 1 - target @ u
    assert deficit <= 4 * epsilon + 1e-12
    print("finite witness extraction geometry: PASS")


def verify_spectral_gap_radius() -> None:
    rng = np.random.default_rng(20260817)
    d, n = 3, 7
    raw = rng.normal(size=(n, d))
    basis, _ = np.linalg.qr(raw)
    projection = basis @ basis.T
    gamma = 0.4
    h = projection + (1 - gamma) * (np.eye(n) - projection)
    for spin in product((-1, 1), repeat=n):
        x = np.asarray(spin, dtype=float)
        q_value = float(x @ h @ x / n)
        defect = 1 - q_value
        u_norm2 = float(x @ projection @ x / n)
        assert 1 - u_norm2 <= defect / gamma + 1e-10
    print("spectral-gap radial certificate: PASS")


def main() -> None:
    verify_circle_covering()
    verify_trust_gap_decomposition()
    verify_finite_extraction_geometry()
    verify_spectral_gap_radius()
    print("Boolean active-eigenspace synchronization checks: PASS")


if __name__ == "__main__":
    main()
