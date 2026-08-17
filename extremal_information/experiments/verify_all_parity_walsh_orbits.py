#!/usr/bin/env python3
"""Exact finite regression for the all-parity rooted Walsh orbit theorem."""

from __future__ import annotations

from itertools import product

import numpy as np


def dot(a: int, b: int) -> int:
    return bin(a & b).count("1") & 1


def apply(matrix: np.ndarray, value: int) -> int:
    m = matrix.shape[0]
    vector = np.array([(value >> j) & 1 for j in range(m)], dtype=np.uint8)
    image = (matrix @ vector) & 1
    return sum(int(image[j]) << j for j in range(m))


def orthogonal_group(m: int) -> list[np.ndarray]:
    identity = np.eye(m, dtype=np.uint8)
    answer = []
    for bits in range(1 << (m * m)):
        matrix = np.array(
            [[(bits >> (i * m + j)) & 1 for j in range(m)] for i in range(m)],
            dtype=np.uint8,
        )
        if np.array_equal((matrix.T @ matrix) & 1, identity):
            answer.append(matrix)
    return answer


def state(labels: tuple[int, ...], m: int) -> tuple[object, ...]:
    k = len(labels)
    gram = tuple(dot(labels[i], labels[j]) for i in range(k) for j in range(k))
    relations = []
    rooted = []
    omega = (1 << m) - 1
    for mask in range(1 << k):
        image = 0
        for i, label in enumerate(labels):
            if (mask >> i) & 1:
                image ^= label
        if image == 0:
            relations.append(mask)
        if image == omega:
            rooted.append(mask)
    return gram, tuple(relations), tuple(rooted)


def canonical_orbit(labels: tuple[int, ...], group: list[np.ndarray]) -> tuple[int, ...]:
    return min(tuple(apply(matrix, label) for label in labels) for matrix in group)


def check_orbit_state(m: int, max_k: int) -> tuple[int, int]:
    group = orthogonal_group(m)
    checked = 0
    classes = 0
    for k in range(1, max_k + 1):
        seen: dict[tuple[object, ...], tuple[int, ...]] = {}
        for labels in product(range(1 << m), repeat=k):
            invariant = state(labels, m)
            orbit = canonical_orbit(labels, group)
            if invariant in seen:
                assert seen[invariant] == orbit
            else:
                seen[invariant] = orbit
            checked += 1
        classes += len(seen)
    return checked, classes


def check_even_group_parameterization_m4() -> int:
    """Construct all T_(S,t,c) after a fixed even-dimensional splitting."""

    m = 4
    # omega=1111, e=1000, and W=<0011,0110>.  In the ordered W basis the
    # restricted form is the order-two symplectic matrix.
    omega = 0b1111
    e = 0b1000
    w_basis = (0b0011, 0b0110)
    assert dot(e, omega) == dot(e, e) == 1
    assert all(dot(e, w) == dot(omega, w) == 0 for w in w_basis)
    assert dot(w_basis[0], w_basis[1]) == 1

    # Enumerate Sp(2,2) directly against J; the ambient dot-product
    # orthogonal group is a different object in characteristic two.
    j = np.array([[0, 1], [1, 0]], dtype=np.uint8)
    sp2 = []
    for bits in range(1 << 4):
        matrix = np.array(
            [[(bits >> (2 * i + h)) & 1 for h in range(2)] for i in range(2)],
            dtype=np.uint8,
        )
        if np.array_equal((matrix.T @ j @ matrix) & 1, j):
            sp2.append(matrix)
    assert len(sp2) == 6

    def w_from_coords(coords: int) -> int:
        value = 0
        for i, basis in enumerate(w_basis):
            if (coords >> i) & 1:
                value ^= basis
        return value

    maps = set()
    source_basis = (e, omega, w_basis[0], w_basis[1])
    # This source basis is invertible.  A linear map is identified by its
    # images on it; compare the resulting standard-coordinate matrices.
    basis_matrix = np.array(
        [[(source_basis[j0] >> i) & 1 for j0 in range(m)] for i in range(m)],
        dtype=np.uint8,
    )
    inverse = np.array(
        [[int(x) for x in row] for row in basis_matrix], dtype=np.uint8
    )
    # Tiny exact Gauss-Jordan inversion over F_2.
    augmented = np.concatenate((inverse, np.eye(m, dtype=np.uint8)), axis=1)
    for col in range(m):
        pivot = next(row for row in range(col, m) if augmented[row, col])
        augmented[[col, pivot]] = augmented[[pivot, col]]
        for row in range(m):
            if row != col and augmented[row, col]:
                augmented[row] ^= augmented[col]
    basis_inverse = augmented[:, m:]

    for symplectic in sp2:
        for t_coords in range(4):
            t = w_from_coords(t_coords)
            for c in range(2):
                images = [e ^ t ^ (omega if c else 0), omega]
                for w_index, w in enumerate(w_basis):
                    sw_coords = sum(
                        int(symplectic[row, w_index]) << row for row in range(2)
                    )
                    sw = w_from_coords(sw_coords)
                    image = sw ^ (omega if dot(t, sw) else 0)
                    images.append(image)
                image_matrix = np.array(
                    [[(images[j0] >> i) & 1 for j0 in range(m)] for i in range(m)],
                    dtype=np.uint8,
                )
                standard = (image_matrix @ basis_inverse) & 1
                assert np.array_equal((standard.T @ standard) & 1, np.eye(m, dtype=np.uint8))
                maps.add(tuple(int(x) for x in standard.reshape(-1)))
    full = orthogonal_group(m)
    assert len(maps) == len(full) == 48
    assert maps == {tuple(int(x) for x in matrix.reshape(-1)) for matrix in full}
    return len(maps)


def sylvester(q: int) -> np.ndarray:
    return np.array([[(-1) ** dot(a, u) for u in range(q)] for a in range(q)], dtype=np.int64)


def modulation(m: int, a: int) -> np.ndarray:
    q = 1 << m
    return np.array([(-1) ** dot(a, v) for u, v in product(range(q), repeat=2)], dtype=np.int64)


def check_even_rooted_collision() -> int:
    m = 4
    q = 1 << m
    n = q * q
    omega = (1 << m) - 1
    other = 0b0011
    assert state((omega,), m)[:2] == state((other,), m)[:2]
    assert state((omega,), m)[2] == (1,)
    assert state((other,), m)[2] == ()

    r = sylvester(q)
    w = np.kron(r, r)
    identity = np.eye(n, dtype=np.int64)
    s0 = np.array(
        [(-1) ** dot(u, v) for u, v in product(range(q), repeat=2)],
        dtype=np.int64,
    )
    s_omega = s0 * modulation(m, omega)
    s_other = s0 * modulation(m, other)
    y_num = w @ s_omega
    assert np.all(y_num % q == 0)
    y = y_num // q
    assert set(y.tolist()) == {-1, 1}

    # The crossed correlation/Rayleigh coordinate is exactly zero.
    h = (s0[:, None] * w) * s0[None, :]
    cross = s_omega * s_other
    assert int(cross @ h @ cross) == 0
    assert np.array_equal(h @ h, n * identity)
    assert np.array_equal(
        (2 * q * identity - h) @ (2 * q * identity + h),
        3 * q * q * identity,
    )
    return 7


def verify() -> None:
    checks_2, classes_2 = check_orbit_state(2, 4)
    checks_4, classes_4 = check_orbit_state(4, 3)
    group_checks = check_even_group_parameterization_m4()
    rooted_checks = check_even_rooted_collision()
    print(
        "all-parity Walsh orbit checks passed: "
        f"{checks_2 + checks_4 + group_checks + rooted_checks} "
        f"({classes_2 + classes_4} orbit/state classes)"
    )


if __name__ == "__main__":
    verify()
