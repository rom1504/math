"""Exact small-field checks for the centered Fourier reflection theorem.

Enumerates every even +/-1 nonconstant multiplier on F_3 and F_3^2,
and every projective Boolean spin. All arithmetic after listing the
characters is integral: N U has kernel sum_line s_l*(2 if a.x=0 else -1).
This checks normalization and finite examples, not the asymptotic proof.
"""

import itertools

import numpy as np


def verify(dimension):
    points = np.array(
        list(itertools.product(range(3), repeat=dimension)), dtype=np.int64
    )
    order = len(points)
    reps = []
    for frequency in points[1:]:
        first = next(int(value) for value in frequency if value)
        if first == 1:
            reps.append(frequency)
    reps = np.array(reps, dtype=np.int64)
    dots = points.dot(reps.T) % 3
    contributions = np.where(dots == 0, 2, -1)
    codes = {tuple(point): index for index, point in enumerate(points)}
    differences = np.array(
        [[codes[tuple((x - y) % 3)] for y in points] for x in points]
    )
    spins = np.ones((2 ** (order - 1), order), dtype=np.int64)
    for bit in range(order - 1):
        spins[:, bit + 1] = 1 - 2 * ((np.arange(len(spins)) >> bit) & 1)
    numerators = []
    for signs in itertools.product((-1, 1), repeat=len(reps)):
        kernel = contributions.dot(np.array(signs, dtype=np.int64))
        matrix = kernel[differences]
        assert np.array_equal(matrix, matrix.T)
        assert np.all(matrix.sum(axis=1) == 0)
        assert np.array_equal(
            matrix.dot(matrix),
            order * order * np.eye(order, dtype=np.int64)
            - order * np.ones((order, order), dtype=np.int64),
        )
        energies = np.einsum("bi,ij,bj->b", spins, matrix, spins)
        numerators.append(int(np.max(np.abs(energies))))
    expected = [8] if dimension == 1 else [72, 80]
    assert sorted(set(numerators)) == expected
    print(
        "PASS F3^%d: %d spectra, %d projective spins each; "
        "absolute Rayleigh numerators %s / %d"
        % (
            dimension,
            len(numerators),
            len(spins),
            sorted(set(numerators)),
            order * order,
        )
    )


if __name__ == "__main__":
    verify(1)
    verify(2)
