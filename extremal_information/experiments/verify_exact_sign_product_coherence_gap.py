#!/usr/bin/env python3
"""Exact checks for the five-port exact-sign coherence gap.

The analytic/scalable statement is in
``drafts/exact_sign_product_coherence_gap.md``.  This verifier checks the
order-16 integer certificate, the depth-two tensor lift, the complete
order-16 Rayleigh shell, and the projective three-port rigidity boundary.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify_bcx_two_port_holonomy import regular_hadamard  # noqa: E402


def word(text: str) -> np.ndarray:
    assert len(text) == 16 and set(text) <= {"+", "-"}
    return np.asarray([1 if symbol == "+" else -1 for symbol in text], dtype=np.int64)


PORTS5 = np.asarray(
    [
        word("-++-+--++-++-++-"),
        word("+++-+-+++--+-+++"),
        word("----------+-----"),
        word("+-----+--++----+"),
        word("+--+-++--+--+--+"),
    ],
    dtype=np.int64,
)
SELECTOR5 = word("+-----+---+----+")

PORTS3 = np.asarray(
    [
        word("+-----+--+-----+"),
        word("+--+--------+--+"),
        word("+++-+-++++-+-+++"),
    ],
    dtype=np.int64,
)
TRIPLE3 = word("+++++--++--+++++")


def cube(n: int):
    return product((-1, 1), repeat=n)


def selector(values) -> int:
    total = sum(values)
    assert total != 0
    return 1 if total > 0 else -1


def selector_witness(ports: np.ndarray, epsilon) -> np.ndarray:
    eps = np.asarray(epsilon, dtype=np.int64)
    return np.asarray(
        [selector(eps * ports[:, coordinate]) for coordinate in range(ports.shape[1])],
        dtype=np.int64,
    )


def odd_subsets(p: int):
    for size in range(1, p + 1, 2):
        yield from combinations(range(p), size)


def port_product(ports: np.ndarray, subset) -> np.ndarray:
    return np.prod(ports[list(subset)], axis=0)


def projective(vector: np.ndarray) -> tuple[int, ...]:
    positive = tuple(int(value) for value in vector)
    negative = tuple(-value for value in positive)
    return min(positive, negative)


def majority_fourier(p: int, subset) -> Fraction:
    numerator = 0
    subset = tuple(subset)
    for values in cube(p):
        character = np.prod([values[i] for i in subset], dtype=int)
        numerator += selector(values) * int(character)
    return Fraction(numerator, 1 << p)


def check_five_port_seed() -> int:
    r, h, hollow = regular_hadamard(2)
    n = len(h)
    assert (r, n) == (4, 16)
    assert np.array_equal(h, h.T)
    assert set(np.unique(h)) == {-1, 1}
    assert np.array_equal(h @ h, 16 * np.eye(16, dtype=np.int64))
    assert np.array_equal(h @ np.ones(n, dtype=np.int64), 4 * np.ones(n, dtype=np.int64))
    assert int(np.trace(h)) == 0
    assert np.array_equal(hollow, h - np.diag(np.diag(h)))

    products = []
    coefficients = []
    checks = 0
    for subset in odd_subsets(5):
        z = port_product(PORTS5, subset)
        coefficient = majority_fourier(5, subset)
        expected = Fraction(3, 8) if len(subset) in (1, 5) else Fraction(-1, 8)
        assert coefficient == expected
        assert int(z @ h @ z) == 52
        products.append(z)
        coefficients.append(coefficient)
        checks += 3

    # All 16 labelled products are distinct (some are antipodal), so this is
    # a genuine affine rank-four product table rather than repeated columns.
    assert len({tuple(int(value) for value in z) for z in products}) == 16

    x = selector_witness(PORTS5, (1, 1, 1, 1, 1))
    assert np.array_equal(x, SELECTOR5)
    assert int(x @ h @ x) == 32

    # Verify the Fourier reconstruction exactly after clearing denominator.
    reconstructed = sum(
        int(8 * coefficient) * z
        for coefficient, z in zip(coefficients, products)
    )
    assert np.array_equal(reconstructed, 8 * x)

    # Exact joint defect: a^T(G-R)a = 1/2, versus diagonal 3/16.
    zmat = np.column_stack(products)
    g_num = zmat.T @ zmat
    r_num = zmat.T @ h @ zmat
    a_num = np.asarray([int(8 * coefficient) for coefficient in coefficients])
    assert int(a_num @ g_num @ a_num) == 64 * n
    assert int(a_num @ r_num @ a_num) == 64 * 32
    assert all(int(g_num[i, i]) == n for i in range(16))
    assert all(int(r_num[i, i]) == 52 for i in range(16))

    # Deleting the trace-zero diagonal preserves every Boolean energy.
    for y_values in cube(n):
        y = np.asarray(y_values, dtype=np.int64)
        assert int(y @ hollow @ y) == int(y @ h @ y)
    checks += 8

    # The next positive Boolean Rayleigh shell below 64 is exactly 52.
    energies = set()
    for y_values in cube(n):
        y = np.asarray(y_values, dtype=np.int64)
        energies.add(int(y @ h @ y))
    below_top = max(value for value in energies if value < r * n)
    assert below_top == 52
    checks += 1 << n
    return checks


def check_tensor_amplification() -> int:
    r, h, _ = regular_hadamard(2)
    one = np.ones(16, dtype=np.int64)
    h2 = np.kron(h, h)
    ports2 = np.asarray([np.kron(w, one) for w in PORTS5], dtype=np.int64)
    x2 = selector_witness(ports2, (1, 1, 1, 1, 1))
    assert np.array_equal(x2, np.kron(SELECTOR5, one))

    checks = 0
    for subset in odd_subsets(5):
        z = port_product(PORTS5, subset)
        z2 = port_product(ports2, subset)
        assert np.array_equal(z2, np.kron(z, one))
        assert int(z2 @ h2 @ z2) == 52 * r * 16
        # More transparently, the normalized ratio remains 52/64.
        assert Fraction(int(z2 @ h2 @ z2), (r * r) * 256) == Fraction(52, 64)
        checks += 3
    assert Fraction(int(x2 @ h2 @ x2), (r * r) * 256) == Fraction(32, 64)

    # Audit the arbitrary-depth scalar formulas without allocating matrices.
    for depth in range(1, 10):
        n = 16**depth
        scale = 4**depth
        assert Fraction(52 * (4 ** (depth - 1)) * (16 ** (depth - 1)), scale * n) == Fraction(52, 64)
        assert Fraction(32 * (4 ** (depth - 1)) * (16 ** (depth - 1)), scale * n) == Fraction(1, 2)
        assert scale * n == n ** 1.5
        checks += 3
    return checks


def check_three_port_rigidity() -> int:
    r, h, _ = regular_hadamard(2)
    assert np.array_equal(np.prod(PORTS3, axis=0), TRIPLE3)
    poles = list(PORTS3) + [TRIPLE3]
    assert all(np.array_equal(h @ z, r * z) for z in poles)

    expected = [
        (0, -1),
        (1, -1),
        (2, -1),
        (3, 1),
        (3, -1),
        (2, 1),
        (1, 1),
        (0, 1),
    ]
    witnesses = []
    for epsilon, (index, sign) in zip(cube(3), expected):
        witness = selector_witness(PORTS3, epsilon)
        assert np.array_equal(witness, sign * poles[index])
        witnesses.append(witness)
    assert {projective(x) for x in witnesses} == {projective(z) for z in poles}

    # This already implies the equality for all diagonal switchings.  Exhaust
    # their projective representatives to reproduce the original observation.
    pole_array = np.asarray(poles, dtype=np.int64)
    witness_array = np.asarray(witnesses, dtype=np.int64)
    for tail in cube(15):
        switch = np.asarray((1,) + tail, dtype=np.int64)
        switched_poles = np.einsum("bi,ij,bj->b", pole_array * switch, h, pole_array * switch)
        switched_witnesses = np.einsum(
            "bi,ij,bj->b", witness_array * switch, h, witness_array * switch
        )
        assert int(np.max(switched_witnesses)) == int(np.max(switched_poles))

    # At PC.3 tensor depth two, the analogous p=5 projective sets are not
    # equal: they have one common element out of sixteen.
    a, b, c = PORTS3
    c1, c2 = a * b, a * c
    one = np.ones(16, dtype=np.int64)
    base = np.kron(a, a)
    generators = [
        np.kron(c1, one),
        np.kron(c2, one),
        np.kron(one, c1),
        np.kron(one, c2),
    ]
    ports5 = np.asarray([base] + [base * generator for generator in generators])
    product_set = {
        projective(port_product(ports5, subset)) for subset in odd_subsets(5)
    }
    witness_set = {
        projective(selector_witness(ports5, epsilon)) for epsilon in cube(5)
    }
    assert len(product_set) == len(witness_set) == 16
    assert len(product_set & witness_set) == 1
    return (1 << 15) + 32


def main() -> None:
    checks = check_five_port_seed()
    checks += check_tensor_amplification()
    checks += check_three_port_rigidity()
    print(f"exact-sign product-coherence checks passed: {checks}")


if __name__ == "__main__":
    main()
