#!/usr/bin/env python3
"""Exact checks for sign-switching visibility synchronization."""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product

import numpy as np
import sympy as sp


def parity(x: int) -> int:
    return bin(int(x)).count("1") & 1


def dot(x: int, y: int) -> int:
    return parity(x & y)


def bent_q(x: int) -> int:
    return parity((x >> 2) & (x & 3))


def polar_q(x: int, y: int) -> int:
    return bent_q(x ^ y) ^ bent_q(x) ^ bent_q(y)


def regular_h16() -> np.ndarray:
    return np.asarray(
        [
            [1 if bent_q(x) ^ dot(x, y) ^ bent_q(y) == 0 else -1
             for y in range(16)]
            for x in range(16)
        ],
        dtype=np.int64,
    )


def switching_translation(t: int) -> np.ndarray:
    assert dot(t, t) == 0
    diagonal = np.asarray(
        [1 if polar_q(x, t) ^ dot(x, t) == 0 else -1 for x in range(16)],
        dtype=np.int64,
    )
    permutation = np.zeros((16, 16), dtype=np.int64)
    for x in range(16):
        permutation[x, x ^ t] = 1
    return permutation @ np.diag(diagonal)


def pc3_poles() -> list[np.ndarray]:
    rows = [
        [1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, -1, -1, -1, -1, 1],
        [1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, 1],
        [1, 1, 1, -1, 1, -1, 1, 1, 1, 1, -1, 1, -1, 1, 1, 1],
    ]
    a, b, c = [np.asarray(row, dtype=np.int64) for row in rows]
    return [a, b, c, a * b * c]


def projectors(p3: np.ndarray, p5: np.ndarray):
    identity = np.eye(16, dtype=np.int64)
    for s3, s5 in product((1, -1), repeat=2):
        # Return the numerator; the projector is numerator/4.
        yield (s3, s5), (identity + s3 * p3) @ (identity + s5 * p5)


def relative_eigenvalues(g: sp.Matrix, d: sp.Matrix) -> list[sp.Rational]:
    # Restrict the generalized pencil to an independent column basis.
    independent = sp.Matrix(g).columnspace()
    basis = sp.Matrix.hstack(*independent)
    g0 = basis.T * g * basis
    d0 = basis.T * d * basis
    values = (g0.inv() * d0).eigenvals()
    result: list[sp.Rational] = []
    for value, multiplicity in values.items():
        result.extend([sp.factor(value)] * multiplicity)
    return sorted(result, key=float)


def check_pc3() -> int:
    h = regular_h16()
    assert np.array_equal(h @ h, 16 * np.eye(16, dtype=np.int64))
    a = h.copy()
    np.fill_diagonal(a, 0)
    assert np.array_equal(a, a.T)
    assert np.all(np.abs(a[np.triu_indices(16, 1)]) == 1)

    lam = sp.symbols("lam")
    charpoly = sp.factor(sp.Matrix(a).charpoly(lam).as_expr())
    target = (lam**2 - 25) * (lam**2 - 17) ** 4 * (lam**2 - 9) ** 3
    assert sp.expand(charpoly - target) == 0

    p3 = switching_translation(3)
    p5 = switching_translation(5)
    identity = np.eye(16, dtype=np.int64)
    assert np.array_equal(p3 @ p3, identity)
    assert np.array_equal(p5 @ p5, identity)
    assert np.array_equal(p3 @ p5, p5 @ p3)
    assert np.array_equal(p3.T @ a @ p3, a)
    assert np.array_equal(p5.T @ a @ p5, a)

    poles = pc3_poles()
    za, zb, zc, ze = poles
    assert np.array_equal(za - zb - zc + ze, np.zeros(16, dtype=np.int64))
    assert np.linalg.matrix_rank(np.stack(poles, axis=1)) == 3

    expected3 = [(-1, 0), (1, 3), (-1, 2), (1, 1)]
    expected5 = [(-1, 2), (-1, 3), (-1, 0), (-1, 1)]
    for generator, expected in ((p3, expected3), (p5, expected5)):
        for index, (sign, target_index) in enumerate(expected):
            assert np.array_equal(generator @ poles[index], sign * poles[target_index])

    weights: dict[str, dict[tuple[int, int], Fraction]] = {}
    for name, z in (("a", za), ("b", zb)):
        weights[name] = {}
        for signature, numerator in projectors(p3, p5):
            # ||Pi z/sqrt(16)||^2 = z^T numerator z /(4*16).
            weights[name][signature] = Fraction(int(z @ numerator @ z), 64)
    assert weights["a"] == {
        (1, 1): 0, (1, -1): 0,
        (-1, 1): Fraction(1, 2), (-1, -1): Fraction(1, 2),
    }
    assert weights["b"] == {
        (1, 1): 0, (1, -1): Fraction(1, 2),
        (-1, 1): Fraction(1, 2), (-1, -1): 0,
    }

    for z in poles:
        assert Fraction(int(z @ a @ z), 5 * 16) == Fraction(4, 5)

    zmat = sp.Matrix(np.stack(poles, axis=1))
    amat = sp.Matrix(a)
    g = zmat.T * zmat / 16
    rayleigh = zmat.T * amat * zmat / 80
    defect = g - rayleigh
    assert relative_eigenvalues(g, defect) == [0, Fraction(2, 5), Fraction(2, 5)]
    assert sp.Matrix(Fraction(2, 5) * g - defect).is_positive_semidefinite

    # Exact check of the Frobenius mismatch constant after deterministic
    # edge flips.  Signed conjugacy is essential in the count.
    for flip_count in range(1, 6):
        changed = a.copy()
        for i, j in list(combinations(range(16), 2))[:flip_count]:
            changed[i, j] *= -1
            changed[j, i] *= -1
        for generator in (p3, p5):
            difference = changed - generator.T @ changed @ generator
            mismatch = sum(
                difference[i, j] != 0 for i, j in combinations(range(16), 2)
            )
            assert int(np.sum(difference * difference)) == 8 * mismatch

    return 38


def walsh_character(s: int, n: int = 8) -> np.ndarray:
    return np.asarray([1 if dot(s, x) == 0 else -1 for x in range(n)], dtype=np.int64)


def check_non_hadamard_cayley() -> int:
    values = np.asarray([0, -1, -1, -1, -1, -1, 1, 1], dtype=np.int64)
    a = np.asarray([[values[x ^ y] for y in range(8)] for x in range(8)], dtype=np.int64)
    assert np.array_equal(a, a.T)
    assert np.all(np.diag(a) == 0)
    assert np.all(np.abs(a[np.triu_indices(8, 1)]) == 1)

    eigenvalues = []
    characters = []
    for s in range(8):
        z = walsh_character(s)
        eigenvalue = int(sum(values[t] * (1 if dot(s, t) == 0 else -1) for t in range(8)))
        assert np.array_equal(a @ z, eigenvalue * z)
        eigenvalues.append(eigenvalue)
        characters.append(z)
    assert eigenvalues == [-3, 1, -3, 1, -3, 1, 5, 1]
    assert not np.array_equal(a @ a, 25 * np.eye(8, dtype=np.int64))

    for t in range(8):
        permutation = np.zeros((8, 8), dtype=np.int64)
        for x in range(8):
            permutation[x, x ^ t] = 1
        assert np.array_equal(permutation.T @ a @ permutation, a)
        for s, z in enumerate(characters):
            sign = 1 if dot(s, t) == 0 else -1
            assert np.array_equal(permutation @ z, sign * z)

    selected = [1, 3, 5, 6, 7]
    zmat = np.stack([characters[s] for s in selected], axis=1)
    g = sp.Matrix(zmat).T * sp.Matrix(zmat) / 8
    rayleigh = sp.Matrix(zmat).T * sp.Matrix(a) * sp.Matrix(zmat) / 40
    defect = g - rayleigh
    assert g == sp.eye(5)
    assert list(defect.diagonal()) == [Fraction(4, 5)] * 3 + [0, Fraction(4, 5)]
    assert relative_eigenvalues(g, defect) == [0] + [Fraction(4, 5)] * 4
    return 32


def check_tensor_visibility() -> int:
    """The PC.3 orbit--character visibility has an exact tensor law."""
    seed_rows = (
        (Fraction(0), Fraction(1, 2), Fraction(1, 2)),
        (Fraction(1, 2), Fraction(1, 2), Fraction(0)),
    )
    checks = 0
    for depth in range(1, 8):
        visibility = Fraction(1)
        for character_word in product(range(3), repeat=depth):
            best = Fraction(0)
            for orbit_word in product(range(2), repeat=depth):
                weight = Fraction(1)
                for orbit, character in zip(orbit_word, character_word):
                    weight *= seed_rows[orbit][character]
                best = max(best, weight)
            assert best == Fraction(1, 2**depth)
            visibility = min(visibility, best)
            checks += 1
        assert visibility == Fraction(1, 2**depth)
        assert visibility**4 == Fraction(1, 16**depth)
        checks += 2
    return checks


def main() -> None:
    checks = check_pc3() + check_non_hadamard_cayley() + check_tensor_visibility()
    print(f"sign-switching visibility checks passed: {checks}")


if __name__ == "__main__":
    main()
