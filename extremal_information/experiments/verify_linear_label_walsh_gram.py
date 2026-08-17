#!/usr/bin/env python3
"""Exact verifier for the linear-label Walsh Gram obstruction."""

from __future__ import annotations

from itertools import product

import numpy as np


def parity_dot(a: int, b: int) -> int:
    return bin(a & b).count("1") & 1


def sylvester(q: int) -> np.ndarray:
    return np.array(
        [[(-1) ** parity_dot(a, u) for u in range(q)] for a in range(q)],
        dtype=np.int64,
    )


def modulation(m: int, a: int) -> np.ndarray:
    q = 1 << m
    return np.array(
        [(-1) ** parity_dot(a, v) for u, v in product(range(q), repeat=2)],
        dtype=np.int64,
    )


def witness(m: int) -> np.ndarray:
    q = 1 << m
    values = []
    for u, v in product(range(q), repeat=2):
        ub = [(u >> j) & 1 for j in range(m)]
        vb = [(v >> j) & 1 for j in range(m)]
        phase = (
            ub[0] * ub[1]
            + vb[0] * vb[1]
            + ub[2] * vb[2]
            + vb[0] + vb[1] + vb[2]
        ) & 1
        for j in range(3, m):
            phase ^= ub[j] * vb[j]
        values.append((-1) ** phase)
    return np.array(values, dtype=np.int64)


def path3() -> np.ndarray:
    return np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]], dtype=np.int64)


def relation_state(labels: tuple[int, ...], m: int) -> tuple[object, ...]:
    k = len(labels)
    gram = tuple(parity_dot(labels[i], labels[j]) for i in range(k) for j in range(k))
    relations = []
    rooted = []
    omega = (1 << m) - 1
    for mask in range(1 << k):
        value = 0
        for i, label in enumerate(labels):
            if (mask >> i) & 1:
                value ^= label
        if value == 0:
            relations.append(mask)
        if value == omega:
            rooted.append(mask)
    return gram, tuple(relations), tuple(rooted)


def orthogonal_group(m: int) -> list[np.ndarray]:
    identity = np.eye(m, dtype=np.uint8)
    group = []
    for mask in range(1 << (m * m)):
        matrix = np.array(
            [[(mask >> (i*m + j)) & 1 for j in range(m)] for i in range(m)],
            dtype=np.uint8,
        )
        if np.array_equal((matrix.T @ matrix) & 1, identity):
            group.append(matrix)
    return group


def apply_binary_matrix(matrix: np.ndarray, value: int) -> int:
    m = matrix.shape[0]
    vector = np.array([(value >> j) & 1 for j in range(m)], dtype=np.uint8)
    image = (matrix @ vector) & 1
    return sum(int(image[j]) << j for j in range(m))


def verify_orbit_classifier_m3() -> int:
    m = 3
    group = orthogonal_group(m)
    assert len(group) > 1
    total_states = 0
    for k in (1, 2, 3):
        state_to_orbits: dict[tuple[object, ...], set[tuple[int, ...]]] = {}
        for labels in product(range(1 << m), repeat=k):
            orbit = min(
                tuple(apply_binary_matrix(matrix, label) for label in labels)
                for matrix in group
            )
            state_to_orbits.setdefault(relation_state(labels, m), set()).add(orbit)
        assert all(len(orbits) == 1 for orbits in state_to_orbits.values())
        total_states += len(state_to_orbits)
    return total_states


def verify() -> None:
    checks = 0
    for m in (3, 4):
        q = 1 << m
        n = q * q
        r = sylvester(q)
        w = np.kron(r, r)
        identity = np.eye(n, dtype=np.int64)
        assert np.array_equal(w @ w, n * identity)
        checks += 1

        a = 0b111
        b = 0b100
        minus_labels = (a, a, a)
        plus_labels = (a, b, a)
        gram_minus = tuple(
            parity_dot(minus_labels[i], minus_labels[j])
            for i in range(3) for j in range(3)
        )
        gram_plus = tuple(
            parity_dot(plus_labels[i], plus_labels[j])
            for i in range(3) for j in range(3)
        )
        assert gram_minus == gram_plus == (1,) * 9
        checks += 1

        da = modulation(m, a)
        db = modulation(m, b)
        ca = (da[:, None] * w) * da[None, :]
        cb = (db[:, None] * w) * db[None, :]
        x = witness(m)

        # The two modulated factor products are respectively self-dual and
        # anti-self-dual.  These identities also prove the child eigenclaims.
        dax = da * x
        dbx = db * x
        assert np.array_equal(w @ dax, q * dax)
        assert np.array_equal(w @ dbx, -q * dbx)
        assert np.array_equal(ca @ x, q * x)
        assert np.array_equal(cb @ x, -q * x)
        checks += 4

        wx = w @ x
        assert np.all(wx % q == 0)
        y = wx // q
        assert set(y.tolist()) == {-1, 1}
        assert np.array_equal(cb @ y, q * y)
        checks += 3

        # The good word (a,b,a) saturates every separate term.
        blocks = (x, y, x)
        children = (ca, cb, ca)
        doubled_energy = sum(int(z @ c @ z) for z, c in zip(blocks, children))
        doubled_energy += 2 * (
            int(blocks[0] @ w @ blocks[1])
            + int(blocks[1] @ w @ blocks[2])
        )
        assert doubled_energy == 7 * n * q
        checks += 1

        # The bad word uses an odd label, hence exact anticommutation.  The
        # three block identities imply M^2=(I+A(P3)^2) tensor I.
        assert np.array_equal(ca @ ca, n * identity)
        assert np.array_equal(w @ ca, -(ca @ w))
        adjacency = path3()
        scalar_square = np.eye(3, dtype=np.int64) + adjacency @ adjacency
        if m == 3:
            global_matrix = np.kron(np.eye(3, dtype=np.int64), ca)
            global_matrix += np.kron(adjacency, w)
            target = n * np.kron(scalar_square, identity)
            assert np.array_equal(global_matrix @ global_matrix, target)
            observed = max(abs(np.linalg.eigvalsh(global_matrix))) / q
            assert abs(observed - np.sqrt(3.0)) < 1e-10
            checks += 2
        checks += 2

        gap = (7 - 3 * np.sqrt(3.0)) / 2
        assert gap > 0.9
        checks += 1

        if m == 3:
            # Gram + relation kernel misses the characteristic root.
            omega = (1 << m) - 1
            unit = 1 << (m - 1)
            state_omega = relation_state((omega,), m)
            state_unit = relation_state((unit,), m)
            assert state_omega[:2] == state_unit[:2]
            assert state_omega[2] == (1,)
            assert state_unit[2] == ()

            # Audit the rho=0 off-pole resolvent input behind the 1/6 gap.
            s0 = np.array(
                [(-1) ** parity_dot(u, v) for u, v in product(range(q), repeat=2)],
                dtype=np.int64,
            )
            mathcal_h = (s0[:, None] * w) * s0[None, :]
            s_omega = s0 * modulation(m, omega)
            s_unit = s0 * modulation(m, unit)
            y_omega_num = w @ s_omega
            assert np.all(y_omega_num % q == 0)
            y_omega = y_omega_num // q
            assert set(y_omega.tolist()) == {-1, 1}
            cross = s_omega * s_unit
            assert int(cross @ mathcal_h @ cross) == 0
            resolvent_numerator = 2 * q * np.eye(n, dtype=np.int64) + mathcal_h
            assert np.array_equal(
                (2 * q * np.eye(n, dtype=np.int64) - mathcal_h)
                @ resolvent_numerator,
                3 * q * q * np.eye(n, dtype=np.int64),
            )
            # Use the correctly matched omega child/query for the exact 3/2 value.
            c_omega = ((modulation(m, omega)[:, None] * w)
                       * modulation(m, omega)[None, :])
            matched = int(s_omega @ c_omega @ s_omega) // 2
            matched += int(s_omega @ w @ y_omega)
            assert matched * 2 == 3 * n * q
            assert (4 * n * q) // 3 < matched
            checks += 9

    classifier_states = verify_orbit_classifier_m3()
    assert classifier_states > 0
    checks += classifier_states

    print(
        "linear-label Walsh Gram/rooted-orbit checks passed: "
        f"{checks} ({classifier_states} m=3 rooted orbit states)"
    )


if __name__ == "__main__":
    verify()
