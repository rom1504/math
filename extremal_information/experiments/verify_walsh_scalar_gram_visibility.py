#!/usr/bin/env python3
"""Exact checks for the unrooted Walsh scalar Gram-visibility theorem."""

from __future__ import annotations

import math

import numpy as np


def parity_dot(a: int, b: int) -> int:
    return bin(a & b).count("1") & 1


def relation_kernel(labels: tuple[int, ...]) -> tuple[int, ...]:
    answer = []
    for mask in range(1 << len(labels)):
        value = 0
        for i, label in enumerate(labels):
            if (mask >> i) & 1:
                value ^= label
        if value == 0:
            answer.append(mask)
    return tuple(answer)


def fwht(values: np.ndarray) -> np.ndarray:
    out = values.astype(np.int64, copy=True)
    width = 1
    while width < len(out):
        for start in range(0, len(out), 2 * width):
            left = out[start : start + width].copy()
            right = out[start + width : start + 2 * width].copy()
            out[start : start + width] = left + right
            out[start + width : start + 2 * width] = left - right
        width *= 2
    return out


def transvection(value: int, axis: int) -> int:
    return value ^ (axis if parity_dot(value, axis) else 0)


def apply_transvections(value: int, axes: tuple[int, ...]) -> int:
    for axis in axes:
        value = transvection(value, axis)
    return value


def normalized_walsh(m: int) -> np.ndarray:
    q = 1 << m
    return np.asarray(
        [[(-1.0) ** parity_dot(a, x) for x in range(q)] for a in range(q)]
    ) / math.sqrt(q)


def reduced_child(f: np.ndarray, label: int) -> np.ndarray:
    diagonal = np.asarray(
        [(-1.0) ** parity_dot(label, x) for x in range(len(f))]
    )
    return diagonal[:, None] * f * diagonal[None, :]


def main() -> int:
    m = 5
    omega = (1 << m) - 1
    good = (0b00011, 0b01100, 0b01111)
    bad = (0b00011, 0b00101, 0b00110)

    assert relation_kernel(good) == relation_kernel(bad) == (0, 7)
    assert all(parity_dot(a, a) == 0 for a in good + bad)
    assert all(parity_dot(good[i], good[j]) == 0 for i in range(3) for j in range(3))
    assert all(
        parity_dot(bad[i], bad[j]) == (i != j)
        for i in range(3)
        for j in range(3)
    )
    good_values = {0, good[0], good[1], good[2]}
    bad_values = {0, bad[0], bad[1], bad[2]}
    assert omega not in good_values and omega not in bad_values

    # The explicit m=5 Witt transport is a product of isotropic
    # transvections.  Coordinates are u_0,...,u_4,v_0,...,v_4.
    axes = (0x44, 0x47, 0x125)
    assert all(parity_dot(axis, axis) == 0 for axis in axes)
    ell_good = (good[0] << m, good[1] << m)
    p1, p2 = 0b00011, 0b00101
    targets = (p1 | (p1 << m), p2 | (p2 << m))
    assert tuple(apply_transvections(x, axes) for x in ell_good) == targets

    columns = tuple(apply_transvections(1 << i, axes) for i in range(2 * m))
    assert all(
        parity_dot(columns[i], columns[j]) == int(i == j)
        for i in range(2 * m)
        for j in range(2 * m)
    )

    # Transport the standard self-dual chirp x_0(u,v)=(-1)^(u.v) back
    # through O.  The length-n Walsh transform is computed exactly by FWHT.
    n = 1 << (2 * m)
    q = 1 << m
    witness = np.empty(n, dtype=np.int8)
    for z in range(n):
        oz = apply_transvections(z, axes)
        u, v = oz & (q - 1), oz >> m
        witness[z] = (-1) ** parity_dot(u, v)
    assert np.array_equal(fwht(witness), q * witness)
    for label in good:
        diagonal = np.asarray(
            [(-1) ** parity_dot(label, z >> m) for z in range(n)],
            dtype=np.int8,
        )
        switched = diagonal * witness
        assert np.array_equal(fwht(switched), q * switched)

    # Verify the reduced Weyl product and the complete triangle spectra.
    f = normalized_walsh(m)
    adjacency = np.ones((3, 3)) - np.eye(3)
    expected_norms = (3.0, (1.0 + math.sqrt(17.0)) / 2.0)
    for labels, flux, expected_norm in (
        (good, 0, expected_norms[0]),
        (bad, 1, expected_norms[1]),
    ):
        children = tuple(reduced_child(f, a) for a in labels)
        product = children[0] @ children[1] @ children[2]
        assert np.allclose(product, ((-1.0) ** flux) * f, atol=2e-12)
        assert all(np.allclose(child @ f, f @ child, atol=2e-12) for child in children)

        carrier = np.kron(adjacency, f)
        for i, child in enumerate(children):
            carrier[i * q : (i + 1) * q, i * q : (i + 1) * q] = child
        eigenvalues = np.linalg.eigvalsh(carrier)
        assert abs(float(np.max(np.abs(eigenvalues))) - expected_norm) < 2e-11

    exact_good = 9.0 / 2.0
    bad_ceiling = 3.0 * (1.0 + math.sqrt(17.0)) / 4.0
    gap = exact_good - bad_ceiling
    assert abs(gap - 3.0 * (5.0 - math.sqrt(17.0)) / 4.0) < 1e-14
    assert gap > 0.657

    print(
        "Walsh scalar Gram visibility checks passed: "
        f"good={exact_good:.12f}, bad_ceiling={bad_ceiling:.12f}, gap={gap:.12f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
