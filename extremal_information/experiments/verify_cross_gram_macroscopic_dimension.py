#!/usr/bin/env python3
"""Finite checks for the cross-Gram macroscopic-dimension draft.

The proof of CG.1--CG.2 is algebraic.  This script checks its two sharpness
models and the total-order normalization used in Sections 3--4:

* a matching cube in the correlation elliptope;
* the Boolean tensor-Walsh realization of a linear number of matching bits;
* the exact order-16 two-port cap gap which exposes one such bit locally.
"""

from __future__ import annotations

from itertools import product
from math import sqrt

import numpy as np

from verify_bounded_cap_contextual_metric_compiler import build


V0 = np.asarray(
    (
        -1, -1, -1, 1,
        -1, -1, 1, -1,
        1, -1, 1, 1,
        -1, 1, 1, 1,
    ),
    dtype=np.int64,
)


def verify_matching_psd_cube() -> None:
    """Enumerate a full affine cube saturating the linear-order regime."""

    p = 6
    eta = 0.75
    edges = ((0, 1), (2, 3), (4, 5))
    for signs in product((-1, 1), repeat=len(edges)):
        matrix = np.eye(p)
        for sign, (i, j) in zip(signs, edges):
            matrix[i, j] = matrix[j, i] = sign * eta
        assert np.linalg.eigvalsh(matrix).min() >= -1e-12

    squared_amplitude = len(edges) * eta**2
    assert squared_amplitude <= p / 2
    print(
        "matching PSD cube:",
        f"h={len(edges)}, sum_eta2={squared_amplitude:.6f}, budget={p/2:.6f}",
    )


def gram_rayleigh(
    operator: np.ndarray, vectors: tuple[np.ndarray, ...], eigenvalue: int
) -> tuple[np.ndarray, np.ndarray]:
    n = vectors[0].size
    gram = np.asarray([[float(a @ b) / n for b in vectors] for a in vectors])
    rayleigh = np.asarray(
        [
            [float(a @ operator @ b) / (eigenvalue * n) for b in vectors]
            for a in vectors
        ]
    )
    return gram, rayleigh


def verify_tensor_walsh_cube() -> None:
    """Realize p/2 independently toggleable cross entries by Boolean ports."""

    q, base_n, h_list, _ = build(2)
    h = np.asarray(h_list, dtype=np.int64)
    one = np.ones(base_n, dtype=np.int64)
    assert one @ V0 == 0
    assert np.array_equal(h @ one, q * one)
    assert np.array_equal(h @ V0, q * V0)

    operator = np.kron(h, h)
    eigenvalue = q * q
    words = (
        np.kron(one, one),
        np.kron(one, V0),
        np.kron(V0, one),
        np.kron(V0, V0),
    )
    n = base_n * base_n
    word_matrix = np.stack(words)
    assert np.array_equal(word_matrix @ word_matrix.T, n * np.eye(4, dtype=int))
    for word in words:
        assert set(word.tolist()) <= {-1, 1}
        assert np.all(np.abs(word) == 1)
        assert np.array_equal(operator @ word, eigenvalue * word)

    observed = {}
    for bit0, bit1 in product((0, 1), repeat=2):
        # A zero selects the orthogonal partner; a one repeats the first word.
        ports = (
            words[0], words[0] if bit0 else words[1],
            words[2], words[2] if bit1 else words[3],
        )
        gram, rayleigh = gram_rayleigh(operator, ports, eigenvalue)
        assert np.allclose(gram, rayleigh)
        assert np.allclose(np.diag(gram), 1)
        expected = np.eye(4)
        expected[0, 1] = expected[1, 0] = bit0
        expected[2, 3] = expected[3, 2] = bit1
        assert np.allclose(gram, expected)
        observed[(bit0, bit1)] = gram

    # In +/- coordinates tau=2*bit-1, each sector has g=r=1/2 per bit.
    p = 4
    h_bits = 2
    gram_rayleigh_energy = h_bits * (0.5**2 + 0.5**2)
    assert gram_rayleigh_energy <= 2 * p
    assert not np.array_equal(observed[(0, 0)], observed[(1, 1)])
    print(
        "tensor Walsh cube:",
        f"n={n}, p={p}, h={h_bits}, sum_(g2+r2)={gram_rayleigh_energy:.6f}",
    )


def hollow_energy(matrix: np.ndarray, spin: tuple[int, ...]) -> int:
    n = len(spin)
    return sum(
        int(matrix[i, j]) * spin[i] * spin[j]
        for i in range(n)
        for j in range(i + 1, n)
    )


def port_cap(matrix: np.ndarray, ports: tuple[np.ndarray, ...], r: int) -> int:
    best = 0
    for spin in product((-1, 1), repeat=matrix.shape[0]):
        value = abs(hollow_energy(matrix, spin))
        x = np.asarray(spin, dtype=np.int64)
        value += r * sum(abs(int(port @ x)) for port in ports)
        best = max(best, value)
    return best


def verify_local_exposure_and_total_scale() -> None:
    q, n, h_list, _ = build(2)
    h = np.asarray(h_list, dtype=np.int64)
    one = np.ones(n, dtype=np.int64)
    repeated = port_cap(h, (one, one), q)
    orthogonal = port_cap(h, (one, V0), q)
    assert repeated == 5 * q * n // 2
    assert orthogonal <= (0.5 + sqrt(2)) * q * n + 1e-9
    assert repeated - orthogonal >= (2 - sqrt(2)) * q * n - 1e-9

    parent_n = 256
    ratios = []
    for p in (1, 4, 16, 64, 256):
        total_n = parent_n + p * int(sqrt(parent_n))
        ratios.append((parent_n / total_n) ** 1.5)
    assert all(left > right for left, right in zip(ratios, ratios[1:]))
    assert ratios[1] > 0.7  # p=n^(1/4) in the tensor example.
    assert ratios[-1] < 0.02
    print(
        "local exposure and total scale:",
        f"caps=({repeated},{orthogonal}), gap={repeated-orthogonal},",
        "ratios=" + ",".join(f"{value:.6f}" for value in ratios),
    )


def main() -> None:
    verify_matching_psd_cube()
    verify_tensor_walsh_cube()
    verify_local_exposure_and_total_scale()
    print("cross-Gram macroscopic-dimension checks: PASS")


if __name__ == "__main__":
    main()
