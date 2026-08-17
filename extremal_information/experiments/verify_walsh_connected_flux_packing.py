#!/usr/bin/env python3
"""Exact structural checks for the connected Walsh flux packing draft."""

from __future__ import annotations

from fractions import Fraction
from itertools import product

import sympy as sp


def dot(x: int, y: int) -> int:
    return bin(x & y).count("1") & 1


def span_value(labels: tuple[int, ...], mask: int) -> int:
    value = 0
    for index, label in enumerate(labels):
        if (mask >> index) & 1:
            value ^= label
    return value


def relation_kernel(labels: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(
        mask
        for mask in range(1 << len(labels))
        if span_value(labels, mask) == 0
    )


def expected_kernel(h: int) -> tuple[int, ...]:
    generators = tuple(0b111 << (3 * i) for i in range(h))
    return tuple(
        sorted(
            span_value(generators, mask)
            for mask in range(1 << len(generators))
        )
    )


def flux_labels(h: int, sigma: tuple[int, ...]) -> tuple[int, ...]:
    labels: list[int] = []
    for i, bit in enumerate(sigma):
        base = 4 * i
        u = (1 << base) | (1 << (base + 1))
        if bit == 0:
            v = (1 << (base + 2)) | (1 << (base + 3))
        else:
            v = (1 << base) | (1 << (base + 2))
        labels.extend((u, v, u ^ v))
    return tuple(labels)


def gram(labels: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    return tuple(tuple(dot(a, b) for b in labels) for a in labels)


def verify_state_cube(max_h: int = 5) -> int:
    checks = 0
    for h in range(1, max_h + 1):
        m = 4 * h + 1
        omega = (1 << m) - 1
        kernel = expected_kernel(h)
        seen_fluxes = set()
        for sigma in product((0, 1), repeat=h):
            labels = flux_labels(h, sigma)
            assert all(label >> (m - 1) == 0 for label in labels)
            assert relation_kernel(labels) == kernel
            assert all(
                span_value(labels, mask) != omega
                for mask in range(1 << len(labels))
            )

            matrix = gram(labels)
            assert all(matrix[i][i] == 0 for i in range(3 * h))
            observed = []
            for i in range(h):
                block = tuple(
                    matrix[3 * i + r][3 * i + s]
                    for r in range(3)
                    for s in range(r + 1, 3)
                )
                assert block in ((0, 0, 0), (1, 1, 1))
                observed.append(block[0])
                for j in range(h):
                    if i == j:
                        continue
                    assert all(
                        matrix[3 * i + r][3 * j + s] == 0
                        for r in range(3)
                        for s in range(3)
                    )
            assert tuple(observed) == sigma
            seen_fluxes.add(tuple(observed))
            checks += 1
        assert len(seen_fluxes) == 1 << h
        checks += 1
    return checks


Vertex = tuple[int, int]
Edge = tuple[Vertex, Vertex]


def canonical_edge(u: Vertex, v: Vertex) -> Edge:
    return (u, v) if u < v else (v, u)


def triangle_edges(h: int) -> set[Edge]:
    return {
        canonical_edge((i, r), (i, s))
        for i in range(h)
        for r in range(3)
        for s in range(r + 1, 3)
    }


def path_connectors(h: int) -> set[Edge]:
    if h == 1:
        return set()
    order = tuple((i, r) for r in range(3) for i in range(h))
    return {
        canonical_edge(order[index], order[index + 1])
        for index in range(len(order) - 1)
    }


def dense_connectors(h: int) -> set[Edge]:
    return {
        canonical_edge((i, r), (j, s))
        for i in range(h)
        for j in range(i + 1, h)
        for r in range(3)
        for s in range(3)
    }


def connected(vertices: set[Vertex], edges: set[Edge]) -> bool:
    if not vertices:
        return True
    neighbors = {vertex: set() for vertex in vertices}
    for u, v in edges:
        neighbors[u].add(v)
        neighbors[v].add(u)
    reached = {next(iter(vertices))}
    frontier = list(reached)
    while frontier:
        u = frontier.pop()
        for v in neighbors[u] - reached:
            reached.add(v)
            frontier.append(v)
    return reached == vertices


def verify_query_graphs(max_h: int = 12) -> int:
    checks = 0
    for h in range(1, max_h + 1):
        vertices = {(i, r) for i in range(h) for r in range(3)}
        triangles = triangle_edges(h)
        path = path_connectors(h)
        dense = dense_connectors(h)
        assert path.isdisjoint(triangles)
        assert dense.isdisjoint(triangles)

        path_graph = path | triangles
        assert connected(vertices, path_graph)
        degrees = {vertex: 0 for vertex in vertices}
        for u, v in path_graph:
            degrees[u] += 1
            degrees[v] += 1
        assert max(degrees.values()) <= 4
        assert len(path) == (0 if h == 1 else 3 * h - 1)

        dense_graph = dense | triangles
        assert connected(vertices, dense_graph)
        assert len(dense) == 9 * h * (h - 1) // 2
        assert len(dense_graph) == 3 * h * (3 * h - 1) // 2

        if h >= 2:
            gamma = Fraction(1, 100 * (h - 1))
            assert 0 < gamma <= 1
            for target in range(h):
                onsite = {
                    (i, r): (Fraction(1) if i == target else gamma)
                    for i in range(h)
                    for r in range(3)
                }
                assert set(onsite) == vertices
                assert all(0 < weight <= 1 for weight in onsite.values())
                edge_weights = {edge: Fraction(1) for edge in path}
                for edge in triangles:
                    gadget = edge[0][0]
                    edge_weights[edge] = (
                        Fraction(1) if gadget == target else gamma
                    )
                assert set(edge_weights) == path_graph
                assert all(0 < weight <= 1 for weight in edge_weights.values())
                perturbation = Fraction(9, 2) * (h - 1) * gamma
                assert perturbation == Fraction(9, 200)
                checks += 1
        else:
            checks += 1
    return checks


def verify_local_spectra_and_gap() -> int:
    adjacency = sp.Matrix([[0, 1, 1], [1, 0, 1], [1, 1, 0]])
    good_sector = sp.eye(3) + adjacency
    one_negative = sp.diag(-1, 1, 1) + adjacency
    three_negative = -sp.eye(3) + adjacency

    assert good_sector.eigenvals() == {sp.Integer(3): 1, sp.Integer(0): 2}
    assert one_negative.charpoly().as_expr().factor() == sp.Symbol("lambda") * (
        sp.Symbol("lambda") ** 2 - sp.Symbol("lambda") - 4
    )
    assert three_negative.eigenvals() == {sp.Integer(1): 1, sp.Integer(-2): 2}

    delta = sp.Rational(3, 4) * (5 - sp.sqrt(17))
    epsilon = sp.Rational(9, 200)
    connected_gap = sp.simplify(delta - 2 * epsilon)
    assert connected_gap == delta - sp.Rational(9, 100)
    assert bool(connected_gap > sp.Rational(1, 2))
    assert abs(float(connected_gap) - 0.5676707807867546) < 1e-15
    return 7


def main() -> None:
    checks = (
        verify_state_cube()
        + verify_query_graphs()
        + verify_local_spectra_and_gap()
    )
    print(f"connected Walsh flux packing checks passed: {checks}")


if __name__ == "__main__":
    main()
