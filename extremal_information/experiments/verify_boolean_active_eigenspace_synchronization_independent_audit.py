#!/usr/bin/env python3
"""Independent checks for Boolean active-eigenspace synchronization.

Adds cap-volume constant checks, exact-eigen algebra, and diagnostics for the
missing d>=2 hypothesis in the finite extraction theorem.
"""

from __future__ import annotations

from itertools import product
from math import gamma, pi, sqrt

import numpy as np


def ball_volume(dimension: int) -> float:
    return pi ** (dimension / 2) / gamma(dimension / 2 + 1)


def check_cap_area_constants() -> int:
    checks = 0
    for d in range(2, 41):
        ratio = 2 * ball_volume(d - 1) / (d * ball_volume(d))
        assert ratio <= 1 + 1e-14
        for epsilon in (1e-5, 0.01, 0.2, 0.5):
            sine = sqrt(2 * epsilon - epsilon * epsilon)
            assert sine <= sqrt(2 * epsilon) + 1e-15
            cap_fraction_upper = ratio * sine ** (d - 1)
            assert cap_fraction_upper <= (2 * epsilon) ** ((d - 1) / 2) + 1e-14
            checks += 1
    return checks


def active_example(d: int):
    """Pair-constant top subspace in R^(2d), with exactly 2^d cube points."""

    n = 2 * d
    basis = np.zeros((n, d))
    for j in range(d):
        basis[2 * j : 2 * j + 2, j] = 1 / sqrt(2)
    projection = basis @ basis.T
    spins = np.asarray(list(product((-1.0, 1.0), repeat=n)))
    q_values = np.einsum("bi,ij,bj->b", spins, projection, spins) / n
    support_vectors = spins @ basis / sqrt(n)
    return projection, spins, q_values, support_vectors


def check_gap_decomposition_and_exact_vertices() -> int:
    checks = 0
    for d in (2, 3):
        _, spins, q_values, support_vectors = active_example(d)
        exact = spins[np.abs(q_values - 1) < 1e-12]
        assert len(exact) == 2**d
        rng = np.random.default_rng(170826 + d)
        beta = 0.8
        for _ in range(200):
            direction = rng.normal(size=d)
            direction /= np.linalg.norm(direction)
            values = q_values / 2 + beta * (support_vectors @ direction)
            index = int(np.argmax(values))
            gap = 0.5 + beta - float(values[index])
            defect = 1 - float(q_values[index])
            angular = 1 - float(support_vectors[index] @ direction)
            assert abs(gap - defect / 2 - beta * angular) < 2e-12
            assert defect <= 2 * gap + 2e-12
            assert angular <= gap / beta + 2e-12
            checks += 1
    return checks


def check_extraction_constant() -> int:
    checks = 0
    for epsilon in (1e-8, 1e-4, 0.01, 0.25):
        distance_squared = ((1 + sqrt(2)) * sqrt(epsilon)) ** 2
        radial_deficit = 2 * epsilon
        resulting_deficit = (radial_deficit + distance_squared) / 2
        assert resulting_deficit < 4 * epsilon + 1e-14
        checks += 1
    return checks


def check_d1_extraction_counterexample() -> int:
    # S^0={-1,+1}.  With epsilon<1/4, a library giving support >=1-4epsilon
    # must contain both signs, while (C/sqrt(epsilon))^(d-1)=1 for d=1.
    epsilon = 0.01
    targets = (-1.0, 1.0)
    witnesses = {-1.0: -1.0, 1.0: 1.0}
    assert all(target * witnesses[target] >= 1 - epsilon for target in targets)
    for singleton in (-1.0, 1.0):
        assert min(target * singleton for target in targets) < 1 - 4 * epsilon
    claimed_exponent_bound = 1.0
    assert len(witnesses) == 2 > claimed_exponent_bound
    return 1


def check_spectral_gap_certificate() -> int:
    rng = np.random.default_rng(281726)
    checks = 0
    for d, n, spectral_gap in ((2, 7, 0.2), (3, 9, 0.7), (4, 10, 1.3)):
        raw = rng.normal(size=(n, d))
        basis, _ = np.linalg.qr(raw)
        projection = basis @ basis.T
        # Put a nonconstant spectrum below 1-gamma on V^perp.
        complement = np.eye(n) - projection
        h = projection + (1 - spectral_gap) * complement
        for spin_word in product((-1, 1), repeat=n):
            spin = np.asarray(spin_word, dtype=float)
            q = float(spin @ h @ spin / n)
            radial = 1 - float(spin @ projection @ spin / n)
            assert radial <= (1 - q) / spectral_gap + 1e-10
            checks += 1
    return checks


def check_exact_eigen_bound() -> int:
    checks = 0
    for d in range(2, 51):
        epsilon_floor = 2 ** (-1 - 2 * d / (d - 1))
        assert epsilon_floor >= 1 / 32 - 1e-15
        assert 2**d >= 1
        checks += 1
    return checks


def main() -> None:
    caps = check_cap_area_constants()
    gaps = check_gap_decomposition_and_exact_vertices()
    extraction = check_extraction_constant()
    d1 = check_d1_extraction_counterexample()
    spectral = check_spectral_gap_certificate()
    exact = check_exact_eigen_bound()
    print(
        "Boolean active-eigenspace independent audit: PASS",
        f"cap_constants={caps}",
        f"gap_exact={gaps}",
        f"extraction={extraction}",
        f"d1_diagnostic={d1}",
        f"spectral={spectral}",
        f"exact_floor={exact}",
    )


if __name__ == "__main__":
    main()
