#!/usr/bin/env python3
"""Exact and matrix-free checks for PC.3 diagonal-switching coherence.

The proof is in ``drafts/pc3_diagonal_switching_coherence.md``.  This
script verifies the finite seed identities exactly and stress-tests the
resulting dimension-free inequality at orders 256 and 4096 without ever
forming the order-4096 Hadamard matrix.
"""

from __future__ import annotations

from itertools import combinations, product
from math import pi, sqrt
from pathlib import Path
import sys

import numpy as np
from sympy import Matrix


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify_bcx_two_port_holonomy import regular_hadamard  # noqa: E402
from verify_exact_sign_product_coherence_gap import PORTS3  # noqa: E402


def odd_subsets(p: int):
    for size in range(1, p + 1, 2):
        yield from combinations(range(p), size)


def product_poles(ports: np.ndarray) -> np.ndarray:
    return np.asarray(
        [np.prod(ports[list(subset)], axis=0) for subset in odd_subsets(len(ports))],
        dtype=np.int64,
    )


def selector_witnesses(ports: np.ndarray) -> np.ndarray:
    witnesses = []
    for epsilon in product((-1, 1), repeat=len(ports)):
        field = np.asarray(epsilon, dtype=np.int64) @ ports
        assert np.all(field != 0)
        x = np.where(field > 0, 1, -1).astype(np.int64)
        representative = tuple(int(v) for v in x)
        opposite = tuple(-v for v in representative)
        if opposite < representative:
            x = -x
            representative = opposite
        if representative not in {tuple(int(v) for v in y) for y in witnesses}:
            witnesses.append(x)
    return np.asarray(witnesses, dtype=np.int64)


def pc3_ports(depth: int) -> np.ndarray:
    a, b, c = PORTS3
    one = np.ones(16, dtype=np.int64)
    relative = (a * b, a * c)
    base = a.copy()
    for _ in range(depth - 1):
        base = np.kron(base, a)
    ports = [base]
    for factor in range(depth):
        for generator in relative:
            pieces = [one] * depth
            pieces[factor] = generator
            lifted = pieces[0]
            for piece in pieces[1:]:
                lifted = np.kron(lifted, piece)
            ports.append(base * lifted)
    return np.asarray(ports, dtype=np.int64)


def apply_tensor_h(vectors: np.ndarray, h: np.ndarray, depth: int) -> np.ndarray:
    """Apply H^(tensor depth) to a batch of row vectors."""
    result = np.asarray(vectors, dtype=np.int64).reshape((-1,) + (16,) * depth)
    for axis in range(1, depth + 1):
        result = np.tensordot(result, h, axes=([axis], [1]))
        # tensordot appends the new coordinate; move it back to its old axis.
        result = np.moveaxis(result, -1, axis)
    return result.reshape((len(vectors), -1))


def seed_identities() -> None:
    r, h, _ = regular_hadamard(2)
    poles = product_poles(PORTS3)
    assert len(poles) == 4
    conjugates_num = np.asarray(
        [z[:, None] * h * z[None, :] for z in poles], dtype=np.int64
    )

    # K is the average of diag(z)(H/r)diag(z), so K=L/16.
    ell = np.sum(conjugates_num, axis=0)
    assert np.array_equal(ell, ell.T)
    assert np.array_equal(ell @ ell @ ell, 16**2 * ell)
    assert Matrix(ell.tolist()).rank() == 6
    assert int(np.trace(ell)) == 0
    eigenvalues = np.linalg.eigvalsh(ell.astype(float) / 16)
    assert np.sum(eigenvalues > 0.5) == 3
    assert np.sum(eigenvalues < -0.5) == 3
    assert np.sum(np.abs(eigenvalues) < 1e-9) == 10

    for numerator in conjugates_num:
        # Each normalized conjugate is a symmetric involution.
        assert np.array_equal(numerator @ numerator, r * r * np.eye(16, dtype=np.int64))


def stress_depth(depth: int, trials: int, seed: int) -> None:
    r0, h, _ = regular_hadamard(2)
    ports = pc3_ports(depth)
    poles = product_poles(ports)
    selectors = selector_witnesses(ports)
    n = 16**depth
    r = r0**depth
    assert poles.shape == selectors.shape == (4**depth, n)

    # Product closure and the selector Fourier identity imply that all these
    # vectors are positive top poles before switching.
    assert np.array_equal(apply_tensor_h(poles, h, depth), r * poles)
    assert np.array_equal(apply_tensor_h(selectors, h, depth), r * selectors)

    if depth == 2:
        # A direct dense regression of the Loewner inequality for every
        # projective selector.  The analytic proof works at every depth.
        seed_poles = product_poles(PORTS3)
        ell = sum(z[:, None] * h * z[None, :] for z in seed_poles)
        k2 = np.kron(ell / 16.0, ell / 16.0)
        h2 = np.kron(h, h) / float(r)
        identity = np.eye(n)
        for x in selectors:
            bx = x[:, None] * h2 * x[None, :]
            certificate = identity - 2 * k2 + bx
            assert np.linalg.eigvalsh(certificate)[0] >= -2e-9

    rng = np.random.default_rng(seed)
    worst_slack_num = None
    for trial in range(trials):
        if trial == 0:
            switch = np.ones(n, dtype=np.int64)
        elif trial == 1:
            switch = np.ones(n, dtype=np.int64)
            switch[: n // 3] = -1
        else:
            switch = rng.choice((-1, 1), size=n).astype(np.int64)

        switched_poles = poles * switch
        switched_selectors = selectors * switch
        hp = apply_tensor_h(switched_poles, h, depth)
        hx = apply_tensor_h(switched_selectors, h, depth)
        pole_energy = np.einsum("bi,bi->b", switched_poles, hp)
        selector_energy = np.einsum("bi,bi->b", switched_selectors, hx)

        # Clearing denominators, d_x <= 2 average_z d_z is exactly
        # q*(rn-E_x) <= 2*sum_z(rn-E_z).
        q = len(poles)
        right = 2 * int(np.sum(r * n - pole_energy))
        slacks = right - q * (r * n - selector_energy)
        assert int(np.min(slacks)) >= 0
        least = int(np.min(slacks))
        worst_slack_num = least if worst_slack_num is None else min(worst_slack_num, least)

    print(
        f"depth={depth}, n={n}, ports={len(ports)}, poles={len(poles)}, "
        f"selectors={len(selectors)}, trials={trials}, least cleared slack={worst_slack_num}"
    )


def diffuse_endpoint_sequence() -> None:
    """Compute (DS.35) from the exact three-atom factor law.

    The local relative-generator pair has values/probabilities
    (1,1):1/4, (1,-1):1/2, (-1,1):1/4.  A tensor transform evaluates all
    4^j raw-character correlations without constructing 16^j rows.
    """
    states = np.asarray(((1, 1), (1, -1), (-1, 1)), dtype=np.int64)
    probabilities = np.asarray((0.25, 0.5, 0.25))
    characters = np.stack(
        (
            np.ones(3),
            states[:, 0],
            states[:, 1],
            states[:, 0] * states[:, 1],
        )
    )
    transform = characters * probabilities
    expected = (
        0.75,
        0.5625,
        0.5,
        0.421875,
        0.3779296875,
        0.349365234375,
        0.3255615234375,
        0.306793212890625,
        0.2899589538574219,
    )
    observed = []
    for depth in range(2, 11):
        alpha = np.asarray([1 if t % 2 == 0 else -1 for t in range(depth)])
        beta = np.ones(depth, dtype=np.int64)
        field = np.ones((3,) * depth, dtype=np.int16)
        for axis in range(depth):
            shape = [1] * depth
            shape[axis] = 3
            field += alpha[axis] * states[:, 0].reshape(shape)
            field += beta[axis] * states[:, 1].reshape(shape)
        response = np.where(field > 0, 1.0, -1.0)
        for axis in range(depth):
            response = np.tensordot(transform, response, axes=([1], [axis]))
            response = np.moveaxis(response, 0, axis)
        observed.append(float(np.max(np.abs(response))))
    assert np.max(np.abs(np.asarray(observed) - np.asarray(expected))) < 1e-14
    print("periodic-endpoint max correlations:", ", ".join(f"{x:.8f}" for x in observed))


def field_scale_sequence() -> None:
    """Check the exact finite laws behind Corollary DS.5.

    Counts are integral because each local PC.3 relative state has one of
    the multiplicities 1, 2, 1 among four equiprobable seed rows.
    """
    counts = {1: 1}
    checkpoints = []
    for depth in range(1, 129):
        if (depth - 1) % 2 == 0:
            local = ((0, 3), (2, 1))
        else:
            local = ((-2, 2), (0, 1), (2, 1))
        updated: dict[int, int] = {}
        for value, multiplicity in counts.items():
            for increment, local_multiplicity in local:
                updated[value + increment] = (
                    updated.get(value + increment, 0)
                    + multiplicity * local_multiplicity
                )
        counts = updated
        total = 4**depth
        mean_num = sum(value * multiplicity for value, multiplicity in counts.items())
        second_num = sum(
            value * value * multiplicity for value, multiplicity in counts.items()
        )
        expected_mean = 1 if depth % 2 == 0 else 1.5
        expected_variance = 7 * depth / 4 - (1 if depth % 2 else 0)
        mean = mean_num / total
        variance = second_num / total - mean * mean
        assert abs(mean - expected_mean) < 1e-14
        assert abs(variance - expected_variance) < 1e-12
        if depth in (8, 16, 32, 64, 128):
            absolute_mean = sum(
                abs(value) * multiplicity for value, multiplicity in counts.items()
            ) / total
            checkpoints.append((depth, absolute_mean / sqrt(depth)))

    target = sqrt(7 / (2 * pi))
    assert abs(checkpoints[-1][1] - target) < 0.01
    print(
        "periodic-field ||h||_1/(N sqrt(j)):",
        ", ".join(f"j={depth}:{ratio:.8f}" for depth, ratio in checkpoints),
        f"limit={target:.8f}",
    )

    # Finite-depth constants in the spherical response separation DS.6.
    depth = 128
    total = 4**depth
    absolute_mean = sum(
        abs(value) * multiplicity for value, multiplicity in counts.items()
    ) / total
    second_moment = sum(
        value * value * multiplicity for value, multiplicity in counts.items()
    ) / total
    rho = absolute_mean / sqrt(second_moment)
    kappa = 0.5
    lam = 0.1
    n = 16**depth
    multiplicity = int(lam * sqrt(n / depth))
    b = multiplicity * sqrt(second_moment) / sqrt(n)
    response_gap_bound = b * (rho - sqrt(1 - rho * rho)) - b * b * rho * rho / (
        2 * kappa
    )
    assert response_gap_bound > 0.014
    print(
        f"labelled spherical response bound at j={depth}: "
        f"rho={rho:.8f}, b={b:.8f}, gap={response_gap_bound:.8f}"
    )


def main() -> None:
    seed_identities()
    stress_depth(depth=2, trials=96, seed=2026081701)
    stress_depth(depth=3, trials=24, seed=2026081702)
    diffuse_endpoint_sequence()
    field_scale_sequence()
    print("PC.3 diagonal-switching coherence checks: PASS")


if __name__ == "__main__":
    main()
