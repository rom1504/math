"""Exact finite checks of the affine Paley ground-state radial law."""

import numpy as np


def verify(prime):
    assert prime % 4 == 1
    character = np.array(
        [0]
        + [1 if pow(a, (prime - 1) // 2, prime) == 1 else -1
           for a in range(1, prime)],
        dtype=np.int64,
    )
    indices = np.arange(prime)
    matrix = character[(indices[:, None] - indices[None, :]) % prime]
    spins = np.ones((2 ** (prime - 1), prime), dtype=np.int64)
    for bit in range(prime - 1):
        spins[:, bit + 1] = 1 - 2 * ((np.arange(len(spins)) >> bit) & 1)
    energies = np.einsum("bi,ij,bj->b", spins, matrix, spins)
    assert int(energies.max()) == -int(energies.min())
    witness = spins[int(np.argmax(energies))]
    ground_energy = int(energies.max())
    covariance_sum = np.zeros((prime, prime), dtype=np.int64)
    orientation_sum = 0
    for scale in range(1, prime):
        for translate in range(prime):
            spin = witness[(scale * indices + translate) % prime]
            orientation = character[scale]
            assert orientation * int(spin.dot(matrix).dot(spin)) == ground_energy
            covariance_sum += orientation * np.outer(spin, spin)
            orientation_sum += orientation
    assert orientation_sum == 0
    assert np.array_equal(covariance_sum, ground_energy * matrix)
    print(
        "PASS prime=%d: Q(core)=%d; exact radial covariance numerator "
        "%d A over denominator %d, all %d affine ground states"
        % (prime, ground_energy // 2, ground_energy,
           prime * (prime - 1), prime * (prime - 1))
    )


if __name__ == "__main__":
    for p in (5, 13, 17):
        verify(p)
