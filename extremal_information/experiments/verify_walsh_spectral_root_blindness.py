#!/usr/bin/env python3
"""Exact checks for the Walsh spectral root-blindness theorem."""

from __future__ import annotations

from itertools import product
import random

import numpy as np
import sympy as sp


def dot(a: int, b: int) -> int:
    return bin(a & b).count("1") & 1


def walsh(m: int) -> np.ndarray:
    q = 1 << m
    return np.array(
        [[(-1) ** dot(x, y) for y in range(q)] for x in range(q)],
        dtype=object,
    )


def child(r: np.ndarray, a: int) -> np.ndarray:
    signs = np.array([(-1) ** dot(a, x) for x in range(r.shape[0])], dtype=object)
    return (signs[:, None] * r) * signs[None, :]


def relation_state(labels: tuple[int, ...]) -> tuple[object, ...]:
    k = len(labels)
    gram = tuple(dot(labels[i], labels[j]) for i in range(k) for j in range(k))
    relations = []
    for mask in range(1 << k):
        value = 0
        for i, label in enumerate(labels):
            if (mask >> i) & 1:
                value ^= label
        if value == 0:
            relations.append(mask)
    return gram, tuple(relations)


def block_carrier(m: int, labels: tuple[int, ...], weights: np.ndarray) -> sp.Matrix:
    q = 1 << m
    r = walsh(m)
    k = len(labels)
    out = np.zeros((k * q, k * q), dtype=object)
    for i, a in enumerate(labels):
        out[i*q:(i+1)*q, i*q:(i+1)*q] = child(r, a)
    for i in range(k):
        for j in range(k):
            if i != j and weights[i, j]:
                out[i*q:(i+1)*q, j*q:(j+1)*q] = int(weights[i, j]) * r
    return sp.Matrix(out.tolist())


def verify_word_trace() -> int:
    checks = 0
    for m, max_length in ((1, 6), (2, 5), (3, 4)):
        q = 1 << m
        r = walsh(m)
        children = [child(r, a) for a in range(q)]
        for length in range(1, max_length + 1):
            for labels in product(range(q), repeat=length):
                matrix = np.eye(q, dtype=object)
                total = 0
                for a in labels:
                    matrix = matrix @ children[a]
                    total ^= a
                # Unnormalised J_a has one factor sqrt(q) per letter.
                trace = int(sum(matrix[i, i] for i in range(q)))
                scale = q ** (length // 2)
                if length % 2 == 0 and total == 0:
                    assert abs(trace) == q * scale
                else:
                    assert trace == 0
                checks += 1
    return checks


def verify_characteristic_polynomials() -> int:
    checks = 0
    cases = [
        (3, (0b111, 0), (0b100, 0)),
        (4, (0b1111, 0), (0b0011, 0)),
    ]
    for m, labels_a, labels_b in cases:
        assert relation_state(labels_a) == relation_state(labels_b)
        for edge_weight in (1, 2, -1):
            weights = np.array([[0, edge_weight], [edge_weight, 0]], dtype=int)
            ka = block_carrier(m, labels_a, weights)
            kb = block_carrier(m, labels_b, weights)
            assert ka.charpoly().all_coeffs() == kb.charpoly().all_coeffs()
            checks += 1

    # Randomly locate further rooted-fibre collisions at m=3 and test two
    # weighted graph extensions exactly.
    rng = random.Random(20260817)
    m = 3
    omega = (1 << m) - 1
    by_state: dict[tuple[object, ...], list[tuple[int, int]]] = {}
    for labels in product(range(1 << m), repeat=2):
        by_state.setdefault(relation_state(labels), []).append(labels)
    tested = 0
    for group in by_state.values():
        roots = []
        for labels in group:
            fibre = tuple(
                mask
                for mask in range(4)
                if ((labels[0] if mask & 1 else 0) ^ (labels[1] if mask & 2 else 0))
                == omega
            )
            roots.append((fibre, labels))
        if len({fibre for fibre, _ in roots}) < 2:
            continue
        first = roots[0]
        second = next(item for item in roots if item[0] != first[0])
        edge_weight = rng.choice((1, 2, -2))
        weights = np.array([[0, edge_weight], [edge_weight, 0]], dtype=int)
        ka = block_carrier(m, first[1], weights)
        kb = block_carrier(m, second[1], weights)
        assert ka.charpoly().all_coeffs() == kb.charpoly().all_coeffs()
        tested += 1
        if tested == 6:
            break
    assert tested > 0
    return checks + tested


def main() -> None:
    checks = verify_word_trace() + verify_characteristic_polynomials()
    print(f"Walsh spectral root-blindness checks passed: {checks}")


if __name__ == "__main__":
    main()
