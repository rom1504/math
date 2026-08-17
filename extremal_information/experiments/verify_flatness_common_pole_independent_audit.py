#!/usr/bin/env python3
"""Independent adversarial checks for the flatness and common-pole drafts.

This script intentionally does not import either primary verifier.  It checks
the identities on an irregular partition tree (including a zero-energy
branch), checks the constants in the recovery estimate on a family whose
spherical optimum is known exactly, and exhausts the order-four common-pole
benchmark.  It also records two scope counterexamples: neither the pair
(||u||_2^2, ||u||_1) nor the scalar synchronization deficit is a complete
contextual response state.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations_with_replacement, product
from math import sqrt

import numpy as np


TOL = 2e-10


def flatness(u: np.ndarray) -> float:
    """Flatness for a vector normalized to squared norm len(u)."""

    return 1.0 - float(np.linalg.norm(u, 1)) / len(u)


def rms(u: np.ndarray) -> float:
    return float(np.linalg.norm(u)) / sqrt(len(u))


# An irregular tree: leaves have different depths.  Each node is a half-open
# interval, and children exactly partition their parent.
TREE = (
    (0, 12),
    (
        ((0, 3), (((0, 1), ()), ((1, 3), ()))),
        ((3, 8), ()),
        (
            (8, 12),
            (
                ((8, 10), (((8, 9), ()), ((9, 10), ()))),
                ((10, 12), ()),
            ),
        ),
    ),
)


def audit_tree(u: np.ndarray) -> tuple[float, dict[int, float]]:
    """Return the FC.2 expansion and transport mass at each depth."""

    root_rms = rms(u)
    n = len(u)
    expansion = 0.0
    depth_mass: dict[int, float] = {}

    def walk(node, depth: int) -> None:
        nonlocal expansion
        (lo, hi), children = node
        block = u[lo:hi]
        block_rms = rms(block)
        relative_rms = block_rms / root_rms
        weight = (hi - lo) / n * relative_rms
        depth_mass[depth] = depth_mass.get(depth, 0.0) + weight

        if not children:
            if block_rms > 0:
                expansion += weight * flatness(block / block_rms)
            return

        # A zero-energy node has zero transport weight.  Its local allocation
        # is immaterial; setting it to zero is the clean convention needed in
        # the theorem statement (relative child RMS values would be 0/0).
        if block_rms == 0:
            local_allocation = 0.0
        else:
            local_allocation = 0.0
            for (child_lo, child_hi), _ in children:
                child = u[child_lo:child_hi]
                lam = (child_hi - child_lo) / (hi - lo)
                rho = rms(child) / block_rms
                local_allocation += 0.5 * lam * (rho - 1.0) ** 2
        expansion += weight * local_allocation
        for child in children:
            walk(child, depth + 1)

    walk(TREE, 0)
    return expansion, depth_mass


def verify_tree_chain_rule() -> int:
    rng = np.random.default_rng(2026081701)
    checks = 0
    for trial in range(250):
        u = rng.normal(size=12)
        if trial % 5 == 0:
            # Exercise the zero-block convention at the first child.
            u[:3] = 0.0
        u *= sqrt(len(u)) / np.linalg.norm(u)
        expansion, depth_mass = audit_tree(u)
        assert abs(expansion - flatness(u)) <= TOL
        assert all(mass <= 1.0 + TOL for mass in depth_mass.values())
        checks += 1 + len(depth_mass)
    return checks


def verify_recovery_constants() -> int:
    """Check FC.3 on rank-one quadratic-plus-linear exposed landscapes."""

    rng = np.random.default_rng(2026081702)
    checks = 0
    for n in (3, 7, 13, 31):
        for _ in range(100):
            z = rng.normal(size=n)
            z /= np.linalg.norm(z)
            lam = float(rng.uniform(0.2, 3.0))
            kappa = float(rng.uniform(0.0, 2.0))
            matrix = lam * np.outer(z, z)
            field = kappa * lam * sqrt(n) * z

            # z is the common top direction of both terms, so this spherical
            # optimum and the Boolean optimum are available in closed form.
            spherical = lam * n * (0.5 + kappa)
            z_l1 = float(np.linalg.norm(z, 1))
            boolean = 0.5 * lam * z_l1**2 + kappa * lam * sqrt(n) * z_l1
            exposed = sqrt(n) * z
            phi = flatness(exposed)
            rhs = lam * n * (1.0 + kappa) * sqrt(2.0 * phi)
            assert -TOL <= spherical - boolean <= rhs + TOL

            signed = np.where(exposed >= 0.0, 1.0, -1.0)
            distance = np.linalg.norm(exposed - signed)
            assert abs(distance**2 - 2.0 * n * phi) <= TOL
            quad_loss = abs(
                0.5 * exposed @ matrix @ exposed
                - 0.5 * signed @ matrix @ signed
            )
            linear_loss = abs(field @ (exposed - signed))
            assert quad_loss <= lam * n * sqrt(2.0 * phi) + TOL
            assert linear_loss <= kappa * lam * n * sqrt(2.0 * phi) + TOL
            checks += 5
    return checks


def verify_pumpable_flatness() -> int:
    checks = 0
    for delta in (0.03, 0.2, 0.65, 0.95):
        rho_plus = sqrt(1.0 + delta)
        rho_minus = sqrt(1.0 - delta)
        s = (rho_plus + rho_minus) / 2.0
        for depth in range(1, 15):
            amplitudes = np.asarray([1.0])
            for _ in range(depth):
                amplitudes = np.kron(amplitudes, (rho_plus, rho_minus))
            n = 2**depth
            assert abs(float(amplitudes @ amplitudes) - n) <= 5e-9
            assert abs(float(np.linalg.norm(amplitudes, 1)) / n - s**depth) <= TOL
            level_sum = sum(s**j * (1.0 - s) for j in range(depth))
            assert abs(level_sum - (1.0 - s**depth)) <= TOL
            checks += 3
    return checks


H4 = np.asarray(
    [
        [1, 1, 1, -1],
        [1, -1, 1, 1],
        [1, 1, -1, 1],
        [-1, 1, 1, 1],
    ],
    dtype=np.int64,
)
X4 = np.asarray(list(product((-1, 1), repeat=4)), dtype=np.int64)
POLE4 = np.ones(4, dtype=np.int64)


def boolean_response(h: np.ndarray, ports: np.ndarray, multiplicity: int) -> int:
    x = np.asarray(list(product((-1, 1), repeat=len(h))), dtype=np.int64)
    quadratic = np.abs(np.einsum("bi,ij,bj->b", x, h, x)) // 2
    linear = multiplicity * np.sum(np.abs(x @ ports.T), axis=1)
    return int(np.max(quadratic + linear))


def verify_common_pole_exhaustively() -> int:
    # All hypotheses, including top Boolean pole, are checked independently.
    assert np.array_equal(H4 @ H4, 4 * np.eye(4, dtype=np.int64))
    assert int(np.trace(H4)) == 0
    assert np.array_equal(H4 @ POLE4, 2 * POLE4)
    checks = 3

    # Exhaust all one- and two-port multisets at the smallest regular order.
    for p in (1, 2):
        for indices in combinations_with_replacement(range(len(X4)), p):
            ports = X4[list(indices)]
            corr_sum = sum(abs(int(w @ POLE4)) for w in ports)
            for multiplicity in (1, 2, 4):
                boolean = boolean_response(H4, ports, multiplicity)
                spherical_upper = 4 + 4 * multiplicity * p
                # CS.5 after cancellation of c and r is exactly mpn*delta.
                deficit_cost = multiplicity * (4 * p - corr_sum)
                assert 0 <= spherical_upper - boolean <= deficit_cost
                checks += 1
    return checks


def verify_tensor_deficit_and_mass() -> int:
    rng = np.random.default_rng(2026081703)
    checks = 0
    for _ in range(100):
        p1, p2 = (int(rng.integers(1, 5)), int(rng.integers(1, 5)))
        w1 = rng.choice((-1, 1), size=(p1, 4)).astype(np.int64)
        w2 = rng.choice((-1, 1), size=(p2, 4)).astype(np.int64)
        c1 = sum(abs(int(w @ POLE4)) for w in w1)
        c2 = sum(abs(int(w @ POLE4)) for w in w2)
        quality1 = Fraction(c1, 4 * p1)
        quality2 = Fraction(c2, 4 * p2)

        tensor_ports = np.asarray([np.kron(a, b) for a in w1 for b in w2])
        tensor_pole = np.kron(POLE4, POLE4)
        tensor_corr = sum(abs(int(w @ tensor_pole)) for w in tensor_ports)
        tensor_quality = Fraction(tensor_corr, 16 * p1 * p2)
        assert tensor_quality == quality1 * quality2

        deficit = 1 - tensor_quality
        deficit1 = 1 - quality1
        deficit2 = 1 - quality2
        assert deficit == deficit1 + deficit2 - deficit1 * deficit2

        # If repetitions tensorize, the other recovery coordinate c=mp/r
        # tensorizes too.  This is not encoded by the deficit alone.
        m1, m2 = int(rng.integers(1, 5)), int(rng.integers(1, 5))
        mass1 = Fraction(m1 * p1, 2)
        mass2 = Fraction(m2 * p2, 2)
        tensor_mass = Fraction((m1 * m2) * (p1 * p2), 4)
        assert tensor_mass == mass1 * mass2
        checks += 3
    return checks


def verify_scope_counterexamples() -> int:
    checks = 0

    # Same (E,L), different value against the same future linear query.
    u = np.asarray([2.0, 1.0, -0.5, -1.5])
    v = np.asarray([1.0, 2.0, -0.5, -1.5])
    assert float(u @ u) == float(v @ v)
    assert float(np.linalg.norm(u, 1)) == float(np.linalg.norm(v, 1))
    query = np.asarray([1.0, -1.0, 0.0, 0.0])
    assert float(query @ u) != float(query @ v)
    checks += 3

    # Same common-pole deficit, different actual Boolean response.  Thus
    # delta is a compositional certificate coordinate, not a complete state.
    ports_a = X4[[0, 3]]
    ports_b = X4[[1, 1]]
    corr_a = sum(abs(int(w @ POLE4)) for w in ports_a)
    corr_b = sum(abs(int(w @ POLE4)) for w in ports_b)
    assert corr_a == corr_b == 4  # p=2,n=4: both deficits are 1/2.
    response_a = boolean_response(H4, ports_a, 1)
    response_b = boolean_response(H4, ports_b, 1)
    assert response_a == 8 and response_b == 10
    checks += 2
    return checks


def main() -> None:
    checks = 0
    checks += verify_tree_chain_rule()
    checks += verify_recovery_constants()
    checks += verify_pumpable_flatness()
    checks += verify_common_pole_exhaustively()
    checks += verify_tensor_deficit_and_mass()
    checks += verify_scope_counterexamples()
    print(f"independent flatness/common-pole audit passed: {checks} checks")


if __name__ == "__main__":
    main()
