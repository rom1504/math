#!/usr/bin/env python3
"""Finite checks for the parity-pole amplification no-go.

The proof is algebraic; this script audits the exact seed data, a padded
odd-monomial selector presentation, tensor/direct-mixture formulas, and the
stable constant in PA.11.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify_exact_sign_product_coherence_gap import (  # noqa: E402
    PORTS5,
    SELECTOR5,
    odd_subsets,
    port_product,
    regular_hadamard,
    selector_witness,
)


def check_seed_shell() -> int:
    _, h, _ = regular_hadamard(2)
    checks = 0
    for subset in odd_subsets(5):
        z = port_product(PORTS5, subset)
        assert Fraction(int(z @ h @ z), 4 * 16) == Fraction(13, 16)
        checks += 1
    assert Fraction(int(SELECTOR5 @ h @ SELECTOR5), 4 * 16) == Fraction(1, 2)

    # The exact top-monomial subgroup used by PA.1b--PA.1c.
    energies = {}
    vectors = {}
    for mask in range(32):
        z = np.ones(16, dtype=np.int64)
        for index in range(5):
            if mask & (1 << index):
                z *= PORTS5[index]
        energies[mask] = int(z @ h @ z)
        vectors[mask] = z
    top = {mask for mask, energy in energies.items() if energy == 64}
    assert top == {0, 5, 9, 12, 17, 20, 24, 29}
    assert {left ^ right for left in top for right in top} == top
    assert {energy: list(energies.values()).count(energy) for energy in set(energies.values())} == {
        32: 4,
        48: 4,
        52: 16,
        64: 8,
    }
    assert all(vectors[mask][3] == vectors[mask][10] for mask in top)
    assert SELECTOR5[3] == -SELECTOR5[10]
    return checks + 6


def check_odd_monomial_padding() -> int:
    """Add a cancelling pair of a degree-three odd monomial.

    The seven-input majority is unchanged, while the full product remains
    an odd seed monomial as predicted by PA.1.
    """

    _, h, _ = regular_hadamard(2)
    triple = PORTS5[0] * PORTS5[1] * PORTS5[2]
    lifted = np.vstack([PORTS5, triple, -triple])
    x = selector_witness(lifted, (1,) * 7)
    assert np.array_equal(x, SELECTOR5)

    full = np.prod(lifted, axis=0)
    expected = -np.prod(PORTS5, axis=0)
    assert np.array_equal(full, expected)
    assert Fraction(int(full @ h @ full), 4 * 16) == Fraction(13, 16)

    # Exhaust every odd subset parity identity for the base five generators.
    checks = 3
    for choices in product(tuple(odd_subsets(5)), repeat=3):
        exponent = set()
        for subset in choices:
            exponent.symmetric_difference_update(subset)
        assert len(exponent) % 2 == 1
        checks += 1
    return checks


def check_even_voters_cannot_make_odd_output() -> int:
    """Audit the parity premise behind PA.1b for all triples of characters."""

    # Character values on the abstract five-cube; this does not depend on
    # whether the 32 patterns occur among the 16 concrete seed rows.
    patterns = list(product((-1, 1), repeat=5))
    even_subsets = [
        subset
        for size in (0, 2, 4)
        for subset in combinations(range(5), size)
    ]

    def character(values, subset):
        answer = 1
        for index in subset:
            answer *= values[index]
        return answer

    checks = 0
    for subsets in product(even_subsets, repeat=3):
        for values in patterns:
            negated = tuple(-value for value in values)
            votes = [character(values, subset) for subset in subsets]
            negated_votes = [character(negated, subset) for subset in subsets]
            assert votes == negated_votes
            assert sum(votes) != 0  # three voters
            assert (sum(votes) > 0) == (sum(negated_votes) > 0)
            checks += 1
    return checks


def check_parity_packing() -> int:
    """Exhaust the affine-hyperplane count behind PA.1c."""

    checks = 0
    for p in (1, 3, 5, 7, 9):
        odd_masks = [mask for mask in range(1 << p) if bin(mask).count("1") % 2]
        for q in range(1, 1 << p):
            odd_images = sum(bin(q & mask).count("1") % 2 for mask in odd_masks)
            if q == (1 << p) - 1:
                assert odd_images == len(odd_masks)
            else:
                assert 2 * odd_images == len(odd_masks)
            average_defect = Fraction(odd_images, len(odd_masks)) * Fraction(3, 16)
            assert average_defect >= Fraction(3, 32)
            checks += 3
    return checks


def check_tensor_formula() -> int:
    _, h, _ = regular_hadamard(2)
    # Use three exact auxiliary contractions with Boolean Rayleigh ratios
    # 1, 1/2, and -1. The first comes from the positive regular pole; the
    # latter two audit the scalar formula independently of exact-sign origin.
    ratios = [Fraction(1), Fraction(1, 2), Fraction(-1)]
    checks = 0
    for s in ratios:
        d_x = 1 - s / 2
        d_z = 1 - Fraction(13, 16) * s
        assert d_z - Fraction(3, 8) * d_x == Fraction(5, 8) * (1 - s)
        assert d_z >= Fraction(3, 8) * d_x
        checks += 2

    # Physical exact-sign top decoration at depth two.
    one = np.ones(16, dtype=np.int64)
    h2 = np.kron(h, h)
    x2 = np.kron(SELECTOR5, one)
    z2 = np.kron(np.prod(PORTS5, axis=0), one)
    assert Fraction(int(x2 @ h2 @ x2), 16 * 256) == Fraction(1, 2)
    assert Fraction(int(z2 @ h2 @ z2), 16 * 256) == Fraction(13, 16)
    checks += 2

    # A recursive outer odd-majority has a nonzero full-leaf coefficient,
    # and its exposed full product becomes progressively worse.
    from math import comb

    for factors in range(1, 13, 2):
        outer_full_coefficient = Fraction(
            (-1) ** ((factors - 1) // 2)
            * comb(factors - 1, (factors - 1) // 2),
            2 ** (factors - 1),
        )
        inner_full_coefficient = Fraction(3, 8)
        assert outer_full_coefficient * inner_full_coefficient**factors != 0
        exposed_ratio = Fraction(13, 16) ** factors
        assert 1 - exposed_ratio >= Fraction(3, 16)
        checks += 2
    return checks


def check_mixtures_and_stability() -> int:
    checks = 0
    # Rational convex mixtures of auxiliary Rayleigh ratios. Linearity must
    # retain PA.9 exactly.
    grid = [Fraction(k, 8) for k in range(-8, 9)]
    for s1 in grid:
        for s2 in grid:
            for weight in [Fraction(k, 8) for k in range(9)]:
                dx = weight * (1 - s1 / 2) + (1 - weight) * (1 - s2 / 2)
                dz = weight * (1 - Fraction(13, 16) * s1) + (
                    1 - weight
                ) * (1 - Fraction(13, 16) * s2)
                assert dz >= Fraction(3, 8) * dx
                checks += 1

    # Exhaust endpoint errors on a rational grid. The 11/8 bound is the
    # triangle-inequality constant in PA.11.
    for s in grid:
        for eta in [Fraction(k, 32) for k in range(9)]:
            for ex_sign, ez_sign in product((-1, 1), repeat=2):
                ex = ex_sign * eta
                ez = ez_sign * eta
                dx = 1 - (s / 2 + ex)
                dz = 1 - (Fraction(13, 16) * s + ez)
                assert dz >= Fraction(3, 8) * dx - Fraction(11, 8) * eta
                checks += 1
    return checks


def main() -> None:
    checks = check_seed_shell()
    checks += check_odd_monomial_padding()
    checks += check_even_voters_cannot_make_odd_output()
    checks += check_parity_packing()
    checks += check_tensor_formula()
    checks += check_mixtures_and_stability()
    print(f"fixed-seed amplification no-go checks passed: {checks}")


if __name__ == "__main__":
    main()
