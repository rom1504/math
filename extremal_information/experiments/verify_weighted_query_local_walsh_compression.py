#!/usr/bin/env python3
"""Exact finite checks for weighted query-local Walsh compression."""

from __future__ import annotations

from itertools import product

import numpy as np


def dot(a: int, b: int) -> int:
    return bin(a & b).count("1") & 1


def walsh(dimension: int) -> np.ndarray:
    size = 1 << dimension
    return np.asarray(
        [[(-1) ** dot(x, y) for y in range(size)] for x in range(size)],
        dtype=np.int16,
    )


def child(w: np.ndarray, label: int, m: int) -> np.ndarray:
    modulation = label << m
    d = np.asarray(
        [(-1) ** dot(modulation, z) for z in range(len(w))], dtype=np.int16
    )
    return (d[:, None] * w) * d[None, :]


def projective_spins(size: int) -> np.ndarray:
    values = np.arange(1 << (size - 1), dtype=np.uint64)[:, None]
    bits = ((values >> np.arange(size - 1, dtype=np.uint64)) & 1).astype(np.int8)
    return np.concatenate(
        [np.ones((len(values), 1), dtype=np.int8), 1 - 2 * bits], axis=1
    ).astype(np.int16)


def canonical_partitions(k: int) -> tuple[tuple[int, ...], ...]:
    answer = []
    for assignment in product(range(k), repeat=k):
        if assignment[0] != 0:
            continue
        seen = 0
        valid = True
        for value in assignment:
            if value > seen + 1:
                valid = False
                break
            seen = max(seen, value)
        if valid:
            answer.append(tuple(assignment))
    return tuple(answer)


def extrema(
    spins: np.ndarray,
    labels: tuple[int, ...],
    vertex_weights: tuple[int, ...],
    edges: tuple[tuple[int, int, int], ...],
    w: np.ndarray,
    children: tuple[np.ndarray, ...],
) -> tuple[int, int]:
    k = len(labels)
    n = len(w)
    blocks = spins.reshape(len(spins), k, n)
    value = np.zeros(len(spins), dtype=np.int64)
    for i, label in enumerate(labels):
        value += (
            vertex_weights[i]
            * np.einsum(
                "bi,ij,bj->b", blocks[:, i], children[label], blocks[:, i],
                optimize=True,
            )
            // 2
        )
    for i, j, weight in edges:
        value += weight * np.einsum(
            "bi,ij,bj->b", blocks[:, i], w, blocks[:, j], optimize=True
        )
    return int(np.min(value)), int(np.max(value))


def verify() -> None:
    m = 1
    w = walsh(2 * m)
    n = len(w)
    children = tuple(child(w, label, m) for label in range(1 << m))
    k = 4
    spins = projective_spins(k * n)
    vertex_weight_families = ((2, -1, 3, 1), (0, -1, 0, 2))
    graph_families = (
        ((0, 1, 1), (1, 2, -2), (2, 3, 1)),
        ((0, 1, -1), (0, 2, 2), (0, 3, -1), (1, 3, 1)),
        ((0, 1, 1), (0, 2, 1), (0, 3, 1), (1, 2, 1), (1, 3, 1), (2, 3, 1)),
    )
    checks = 0
    for labels in product(range(1 << m), repeat=k):
        for vertex_weights in vertex_weight_families:
            for edges in graph_families:
                full_min, full_max = extrema(
                    spins, labels, vertex_weights, edges, w, children
                )
                for partition in canonical_partitions(k):
                    kept = tuple(
                        edge for edge in edges if partition[edge[0]] == partition[edge[1]]
                    )
                    deleted_mass = sum(
                        abs(edge[2])
                        for edge in edges
                        if partition[edge[0]] != partition[edge[1]]
                    )
                    truncated_min, truncated_max = extrema(
                        spins, labels, vertex_weights, kept, w, children
                    )
                    bound = deleted_mass * n ** 1.5 + 1e-9
                    assert abs(full_max - truncated_max) <= bound
                    assert abs(full_min - truncated_min) <= bound
                    assert abs(
                        max(full_max, -full_min)
                        - max(truncated_max, -truncated_min)
                    ) <= bound
                    checks += 3

    # Pure combinatorial checks for the path/dense regimes.
    for t in range(2, 41):
        for blocks in canonical_partitions(min(t, 7)) if t <= 7 else ():
            sizes = [blocks.count(i) for i in set(blocks)]
            cross = (t * t - sum(s * s for s in sizes)) // 2
            if t <= 7:
                direct = sum(
                    1
                    for i in range(t)
                    for j in range(i + 1, t)
                    if blocks[i] != blocks[j]
                )
                assert cross == direct
                checks += 1
    print(f"weighted query-local Walsh checks passed: {checks}")


if __name__ == "__main__":
    verify()
