#!/usr/bin/env python3
"""Exact verifier for the unrooted Walsh ambient-orbit collapse.

The script constructs explicit orthogonal maps by products of isotropic
transvections, checks the odd/even characteristic-root collisions, and
verifies simultaneous conjugacy of a nontrivial weighted three-block matrix.
"""

from __future__ import annotations

from collections import deque
from itertools import product

import numpy as np


def parity_dot(a: int, b: int) -> int:
    return bin(a & b).count("1") & 1


def bits(value: int, dimension: int) -> np.ndarray:
    return np.asarray([(value >> j) & 1 for j in range(dimension)], dtype=np.uint8)


def integer(vector: np.ndarray) -> int:
    return sum(int(vector[j]) << j for j in range(len(vector)))


def transvection(value: int, vector: int) -> int:
    return value ^ (vector if parity_dot(value, vector) else 0)


def transvection_path(source: int, target: int, dimension: int) -> list[int]:
    """Find isotropic transvections carrying source to target."""

    isotropic = [v for v in range(1, 1 << dimension) if parity_dot(v, v) == 0]
    queue = deque([source])
    parent: dict[int, tuple[int, int] | None] = {source: None}
    while queue:
        value = queue.popleft()
        if value == target:
            break
        for vector in isotropic:
            image = transvection(value, vector)
            if image not in parent:
                parent[image] = (value, vector)
                queue.append(image)
    if target not in parent:
        raise AssertionError((source, target, dimension))
    path: list[int] = []
    value = target
    while parent[value] is not None:
        previous, vector = parent[value]
        path.append(vector)
        value = previous
    return list(reversed(path))


def orthogonal_map(source: int, target: int, dimension: int) -> tuple[np.ndarray, list[int]]:
    path = transvection_path(source, target, dimension)
    matrix = np.eye(dimension, dtype=np.uint8)
    identity = np.eye(dimension, dtype=np.uint8)
    for value in path:
        vector = bits(value, dimension)[:, None]
        step = (identity + vector @ vector.T) & 1
        matrix = (step @ matrix) & 1
    assert integer((matrix @ bits(source, dimension)) & 1) == target
    assert np.array_equal((matrix.T @ matrix) & 1, identity)
    return matrix, path


def apply(matrix: np.ndarray, value: int) -> int:
    return integer((matrix @ bits(value, len(matrix))) & 1)


def orthogonal_group(dimension: int) -> list[np.ndarray]:
    identity = np.eye(dimension, dtype=np.uint8)
    answer: list[np.ndarray] = []
    for mask in range(1 << (dimension * dimension)):
        matrix = np.asarray(
            [
                [(mask >> (i * dimension + j)) & 1 for j in range(dimension)]
                for i in range(dimension)
            ],
            dtype=np.uint8,
        )
        if np.array_equal((matrix.T @ matrix) & 1, identity):
            answer.append(matrix)
    return answer


def walsh(dimension: int) -> np.ndarray:
    size = 1 << dimension
    return np.asarray(
        [[(-1) ** parity_dot(x, y) for y in range(size)] for x in range(size)],
        dtype=np.int16,
    )


def child(matrix: np.ndarray, modulation: int) -> np.ndarray:
    diagonal = np.asarray(
        [(-1) ** parity_dot(modulation, z) for z in range(len(matrix))], dtype=np.int16
    )
    return (diagonal[:, None] * matrix) * diagonal[None, :]


def relation_state(labels: tuple[int, ...], m: int) -> tuple[object, ...]:
    gram = tuple(parity_dot(a, b) for a in labels for b in labels)
    zero: list[int] = []
    root: list[int] = []
    omega = (1 << m) - 1
    for mask in range(1 << len(labels)):
        value = 0
        for i, label in enumerate(labels):
            if (mask >> i) & 1:
                value ^= label
        if value == 0:
            zero.append(mask)
        if value == omega:
            root.append(mask)
    return gram, tuple(zero), tuple(root)


def full_matrix(
    w: np.ndarray,
    children: tuple[np.ndarray, ...],
    vertex_weights: tuple[int, ...],
    edge_weights: np.ndarray,
) -> np.ndarray:
    k = len(children)
    blocks: list[list[np.ndarray]] = []
    for i in range(k):
        row: list[np.ndarray] = []
        for j in range(k):
            if i == j:
                row.append(vertex_weights[i] * children[i])
            else:
                row.append(int(edge_weights[i, j]) * w)
        blocks.append(row)
    return np.block(blocks)


def verify_case(m: int, a: int, b: int) -> int:
    dimension = 2 * m
    omega = (1 << m) - 1
    ambient_omega = omega | (omega << m)
    source = a << m
    target = b << m
    assert parity_dot(source, source) == parity_dot(target, target)
    assert source not in (0, ambient_omega) and target not in (0, ambient_omega)
    orthogonal, path = orthogonal_map(source, target, dimension)
    assert apply(orthogonal, ambient_omega) == ambient_omega

    w = walsh(dimension)
    permutation = np.asarray([apply(orthogonal, z) for z in range(len(w))])
    assert len(set(int(x) for x in permutation)) == len(w)
    assert np.array_equal(w[np.ix_(permutation, permutation)], w)
    ca = child(w, source)
    cb = child(w, target)
    assert np.array_equal(cb[np.ix_(permutation, permutation)], ca)

    k = 3
    labels_a = (a,) * k
    labels_b = (b,) * k
    state_a = relation_state(labels_a, m)
    state_b = relation_state(labels_b, m)
    assert state_a[:2] == state_b[:2]
    assert state_a[2] != state_b[2]

    vertex_weights = (2, -1, 3)
    edge_weights = np.asarray([[0, 1, -2], [1, 0, 4], [-2, 4, 0]], dtype=np.int16)
    ma = full_matrix(w, (ca,) * k, vertex_weights, edge_weights)
    mb = full_matrix(w, (cb,) * k, vertex_weights, edge_weights)
    block_permutation = np.concatenate([i * len(w) + permutation for i in range(k)])
    assert np.array_equal(mb[np.ix_(block_permutation, block_permutation)], ma)

    # Pointwise Boolean energy checks on deterministic, non-symmetric tests.
    rng = np.random.default_rng(1000 + m)
    for _ in range(16):
        x = rng.choice(np.asarray([-1, 1], dtype=np.int16), size=k * len(w))
        px = x[block_permutation]
        assert int(x @ ma @ x) == int(px @ mb @ px)
    return 5 + len(path) + 16


def verify() -> None:
    checks = 0
    # Odd label dimension: omega and e_1 are nonzero anisotropic labels.
    checks += verify_case(m=3, a=0b111, b=0b001)
    # Even label dimension: omega and e_1+e_2 are nonzero isotropic labels.
    checks += verify_case(m=4, a=0b1111, b=0b0011)

    # Exhaust the complete m=2 tuple space through length three.  Equality
    # of (Gram, relation kernel) must coincide with the ambient O(4,2) orbit
    # after embedding a -> (0,a); no characteristic-root field is available.
    m = 2
    group = orthogonal_group(2 * m)
    ambient_omega = ((1 << m) - 1) | (((1 << m) - 1) << m)
    for k in (1, 2, 3):
        state_to_orbits: dict[tuple[object, ...], set[tuple[int, ...]]] = {}
        for labels in product(range(1 << m), repeat=k):
            embedded = tuple(label << m for label in labels)
            span = set()
            for mask in range(1 << k):
                value = 0
                for i in range(k):
                    if (mask >> i) & 1:
                        value ^= embedded[i]
                span.add(value)
            assert ambient_omega not in span
            orbit = min(tuple(apply(matrix, value) for value in embedded) for matrix in group)
            state = relation_state(tuple(labels), m)[:2]
            state_to_orbits.setdefault(state, set()).add(orbit)
        assert all(len(orbits) == 1 for orbits in state_to_orbits.values())
        checks += len(state_to_orbits)

    # Complete tiny Boolean wind tunnel: every m=1, k=3 tuple is separated by
    # the vector of unweighted labelled graph maxima.  This is deliberately a
    # scalar finite observation, independent of the conjugacy theorem above.
    from walsh_semantic_wind_tunnel import m1_signatures

    signatures = m1_signatures(3)
    graph_signatures = {values["graphs"] for values in signatures.values()}
    assert len(signatures) == len(graph_signatures) == 8
    checks += len(signatures)
    print(
        "unrooted Walsh ambient-collapse checks passed: "
        f"{checks} (|O(4,2)|={len(group)})"
    )


if __name__ == "__main__":
    verify()
