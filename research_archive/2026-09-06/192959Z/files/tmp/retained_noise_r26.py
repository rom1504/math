#!/usr/bin/env python3
"""Finite checks for the Wave 26 retained-deficit noise identities.

All expectations are enumerated directly.  The tests deliberately use
ordinary (not projectively quotiented) spin outcomes; quadratic energies and
row-square costs are invariant under the resulting duplicate labels.
"""

from itertools import combinations, product

import numpy as np

from check_response_dual_r16 import A8, A9, qnorm, spins
from verify_compatible_replacement_r12 import A6


def child_ground(a, selected):
    """Return (orientation, spin, Q) for one oriented child ground."""
    block = a[np.ix_(selected, selected)]
    best = None
    for y in spins(len(selected)):
        energy = int(y @ block @ y)
        for orientation in (-1, 1):
            value = orientation * energy
            if best is None or value > best[2]:
                best = (orientation, y.copy(), value)
    return best


def probability(bits, flip_probability):
    negatives = sum(bit == -1 for bit in bits)
    return flip_probability**negatives * (1.0 - flip_probability) ** (
        len(bits) - negatives
    )


def check_child_channel(a, m, flip_probability=0.2, orientation_probability=0.1):
    n = len(a)
    q_parent = qnorm(a)
    theta = 1.0 - 2.0 * flip_probability
    eta = 1.0 - 2.0 * orientation_probability
    retained = eta * theta**2
    rho = m / n
    p2 = m * (m - 1) / (n * (n - 1))
    baseline = (rho**1.5 - p2) * q_parent

    for selected_tuple in combinations(range(n), m):
        selected = np.array(selected_tuple, dtype=int)
        outside = np.array(
            [index for index in range(n) if index not in selected_tuple], dtype=int
        )
        sigma, y, q_child = child_ground(a, selected)

        total_mass = 0.0
        mean_child = 0.0
        mean_parent = 0.0
        mean_loss = 0.0
        mean_deficit = 0.0
        mean_effective = 0.0
        mean_row_square = 0.0

        for noise_tuple in product((-1, 1), repeat=m):
            noise = np.array(noise_tuple, dtype=np.int64)
            noise_mass = probability(noise_tuple, flip_probability)
            for tau in (-1, 1):
                tau_mass = (
                    orientation_probability
                    if tau == -1
                    else 1.0 - orientation_probability
                )
                orientation = sigma * tau
                for outside_tuple in product((-1, 1), repeat=n - m):
                    x = np.empty(n, dtype=np.int64)
                    x[selected] = y * noise
                    x[outside] = np.array(outside_tuple, dtype=np.int64)
                    mass = noise_mass * tau_mass / (2 ** (n - m))

                    child_energy = int(
                        orientation
                        * x[selected]
                        @ a[np.ix_(selected, selected)]
                        @ x[selected]
                    )
                    parent_energy = int(orientation * x @ a @ x)
                    loss = q_child - child_energy
                    deficit = q_parent - parent_energy
                    effective = loss - p2 * deficit - baseline
                    row_square = int(np.sum((a @ x) ** 2))

                    total_mass += mass
                    mean_child += mass * child_energy
                    mean_parent += mass * parent_energy
                    mean_loss += mass * loss
                    mean_deficit += mass * deficit
                    mean_effective += mass * effective
                    mean_row_square += mass * row_square

        row_vector = a[:, selected] @ y
        predicted_row_square = n * (n - 1) + theta**2 * (
            int(row_vector @ row_vector) - m * (n - 1)
        )
        assert abs(total_mass - 1.0) < 1e-11
        assert abs(mean_child - retained * q_child) < 1e-9
        assert abs(mean_parent - retained * q_child) < 1e-9
        assert abs(mean_loss - (1.0 - retained) * q_child) < 1e-9
        assert abs(mean_deficit - (q_parent - retained * q_child)) < 1e-9
        assert abs(
            mean_effective
            - ((1.0 - retained * (1.0 - p2)) * q_child - rho**1.5 * q_parent)
        ) < 1e-9
        assert abs(mean_row_square - predicted_row_square) < 1e-9


def parent_ground(a):
    best = None
    for x in spins(len(a)):
        energy = int(x @ a @ x)
        for orientation in (-1, 1):
            value = orientation * energy
            if best is None or value > best[2]:
                best = (orientation, x.copy(), value)
    return best


def check_parent_noise(a, flip_probability=0.2):
    n = len(a)
    theta = 1.0 - 2.0 * flip_probability
    _, ground, _ = parent_ground(a)
    ground_cost = int(np.sum((a @ ground) ** 2))
    observed = 0.0
    for noise_tuple in product((-1, 1), repeat=n):
        noise = np.array(noise_tuple, dtype=np.int64)
        mass = probability(noise_tuple, flip_probability)
        x = ground * noise
        observed += mass * int(np.sum((a @ x) ** 2))
    predicted = n * (n - 1) + theta**2 * (ground_cost - n * (n - 1))
    assert abs(observed - predicted) < 1e-9


if __name__ == "__main__":
    for matrix in (A6, A8, A9):
        check_child_channel(matrix.astype(np.int64), len(matrix) - 1)
        check_parent_noise(matrix.astype(np.int64))
    print("retained-deficit noise identities: all finite checks passed")
