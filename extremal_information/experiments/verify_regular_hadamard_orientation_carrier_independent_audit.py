#!/usr/bin/env python3
"""Independent finite audit of the regular-Hadamard orientation carrier.

This deliberately adds checks not present in the canonical verifier:

* OC.9 classifies projective-switching orbits for disconnected supports;
* OC.2 has uniform 2^r fibres for non-singleton connected pieces;
* both triangle bridge-product classes have the claimed spectra/caps at n=4.
"""

from __future__ import annotations

import itertools
from collections import defaultdict, deque

import numpy as np


def signs(k: int):
    return itertools.product((-1, 1), repeat=k)


def all_edges(k: int) -> list[tuple[int, int]]:
    return list(itertools.combinations(range(k), 2))


def scalar_matrix(
    sigma: tuple[int, ...],
    active_edges: tuple[tuple[int, int], ...],
    bridge: tuple[int, ...],
) -> np.ndarray:
    matrix = np.diag(np.asarray(sigma, dtype=int))
    for value, (i, j) in zip(bridge, active_edges):
        matrix[i, j] = matrix[j, i] = value
    return matrix


def orbit_key(matrix: np.ndarray) -> tuple[int, ...]:
    k = len(matrix)
    words = []
    for switch_word in signs(k):
        switch = np.asarray(switch_word, dtype=int)
        for antipode in (-1, 1):
            image = antipode * switch[:, None] * matrix * switch[None, :]
            words.append(tuple(int(value) for value in image.ravel()))
    return min(words)


def forest_and_chords(
    k: int, active_edges: tuple[tuple[int, int], ...]
) -> tuple[list[tuple[int, int]], list[tuple[int, int]], list[int]]:
    parent = list(range(k))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    forest = []
    chords = []
    for i, j in active_edges:
        ri, rj = find(i), find(j)
        if ri == rj:
            chords.append((i, j))
        else:
            parent[ri] = rj
            forest.append((i, j))

    components: dict[int, list[int]] = defaultdict(list)
    for i in range(k):
        components[find(i)].append(i)
    roots = sorted(min(component) for component in components.values())
    return forest, chords, roots


def forest_path(
    k: int, forest: list[tuple[int, int]], source: int, target: int
) -> list[tuple[int, int]]:
    adjacency: list[list[int]] = [[] for _ in range(k)]
    for i, j in forest:
        adjacency[i].append(j)
        adjacency[j].append(i)
    queue = deque([source])
    predecessor = {source: -1}
    while queue:
        i = queue.popleft()
        if i == target:
            break
        for j in adjacency[i]:
            if j not in predecessor:
                predecessor[j] = i
                queue.append(j)
    assert target in predecessor
    path = []
    cursor = target
    while predecessor[cursor] != -1:
        previous = predecessor[cursor]
        path.append(tuple(sorted((cursor, previous))))
        cursor = previous
    return path


def canonical_coordinates(
    matrix: np.ndarray, active_edges: tuple[tuple[int, int], ...]
) -> tuple[int, ...]:
    """Compute the invariant coordinates displayed in OC.9."""

    k = len(matrix)
    forest, chords, roots = forest_and_chords(k, active_edges)
    distinguished = roots[0]
    sigma_root = int(matrix[distinguished, distinguished])
    output = [
        int(matrix[v, v]) * sigma_root
        for v in range(k)
        if v != distinguished
    ]
    edge_sign = {
        tuple(sorted((i, j))): int(matrix[i, j]) for i, j in active_edges
    }
    for i, j in chords:
        cycle = forest_path(k, forest, i, j) + [tuple(sorted((i, j)))]
        product_sign = 1
        for edge in cycle:
            product_sign *= edge_sign[edge]
        output.append((sigma_root ** len(cycle)) * product_sign)
    return tuple(output)


def check_canonical_coordinates() -> int:
    checks = 0
    for k in range(1, 5):
        candidate_edges = all_edges(k)
        for support in itertools.product((0, 1), repeat=len(candidate_edges)):
            active_edges = tuple(
                edge for keep, edge in zip(support, candidate_edges) if keep
            )
            coordinate_to_orbit: dict[tuple[int, ...], tuple[int, ...]] = {}
            orbit_to_coordinate: dict[tuple[int, ...], tuple[int, ...]] = {}
            for sigma in signs(k):
                for bridge in signs(len(active_edges)):
                    matrix = scalar_matrix(sigma, active_edges, bridge)
                    coordinate = canonical_coordinates(matrix, active_edges)
                    orbit = orbit_key(matrix)
                    assert coordinate_to_orbit.setdefault(coordinate, orbit) == orbit
                    assert orbit_to_coordinate.setdefault(orbit, coordinate) == coordinate

            _, _, roots = forest_and_chords(k, active_edges)
            expected = 2 ** (len(active_edges) + len(roots) - 1)
            assert len(coordinate_to_orbit) == expected
            assert len(orbit_to_coordinate) == expected
            checks += 1
    return checks


def induced_key(matrix: np.ndarray, vertices: tuple[int, ...]) -> tuple[int, ...]:
    return orbit_key(matrix[np.ix_(vertices, vertices)])


def check_one_gluing_case(
    pieces: tuple[tuple[int, ...], ...],
    internal_edges: tuple[tuple[int, int], ...],
    cross_edges: tuple[tuple[int, int], ...],
) -> int:
    k = sum(len(piece) for piece in pieces)
    active_edges = internal_edges + cross_edges
    whole_to_marginal: dict[tuple[int, ...], tuple[tuple[int, ...], ...]] = {}
    fibres: dict[tuple[tuple[int, ...], ...], set[tuple[int, ...]]] = defaultdict(set)
    for sigma in signs(k):
        for bridge in signs(len(active_edges)):
            matrix = scalar_matrix(sigma, active_edges, bridge)
            whole = orbit_key(matrix)
            marginal = tuple(induced_key(matrix, piece) for piece in pieces)
            assert whole_to_marginal.setdefault(whole, marginal) == marginal
            fibres[marginal].add(whole)

    expected_marginals = 2 ** len(internal_edges)
    assert len(fibres) == expected_marginals
    expected_fibre = 2 ** len(cross_edges)
    assert all(len(fibre) == expected_fibre for fibre in fibres.values())
    return len(fibres)


def check_general_gluing_fibres() -> int:
    checks = 0
    # A two-vertex piece joined to a singleton by one edge.
    checks += check_one_gluing_case(
        ((0, 1), (2,)), ((0, 1),), ((1, 2),)
    )
    # Two non-singleton pieces with two attachment vertices and a cross cycle.
    checks += check_one_gluing_case(
        ((0, 1), (2, 3)), ((0, 1), (2, 3)), ((1, 2), (0, 3))
    )
    # A triangular piece and singleton with two cross edges.
    checks += check_one_gluing_case(
        ((0, 1, 2), (3,)), ((0, 1), (1, 2), (0, 2)), ((0, 3), (2, 3))
    )
    return checks


def exact_cap(matrix: np.ndarray, hadamard: np.ndarray) -> int:
    full = np.kron(matrix, hadamard)
    best = 0
    for spin_word in signs(full.shape[0]):
        spin = np.asarray(spin_word, dtype=int)
        best = max(best, abs(int(spin @ full @ spin)) // 2)
    return best


def check_triangle_classes() -> int:
    hadamard = np.asarray(
        (
            (1, 1, 1, -1),
            (1, -1, 1, 1),
            (1, 1, -1, 1),
            (-1, 1, 1, 1),
        ),
        dtype=int,
    )
    triangle_edges = ((0, 1), (0, 2), (1, 2))
    caps_by_product: dict[int, set[int]] = defaultdict(set)
    spectra_by_product: dict[int, set[tuple[float, ...]]] = defaultdict(set)
    for bridge in signs(3):
        matrix = scalar_matrix((1, 1, 1), triangle_edges, bridge)
        product_sign = bridge[0] * bridge[1] * bridge[2]
        caps_by_product[product_sign].add(exact_cap(matrix, hadamard))
        spectrum = tuple(np.round(np.linalg.eigvalsh(matrix), 9))
        spectra_by_product[product_sign].add(spectrum)
    assert caps_by_product[1] == {36}
    assert caps_by_product[-1] == {20}
    assert spectra_by_product[1] == {(-0.0, 0.0, 3.0)}
    assert spectra_by_product[-1] == {(-1.0, 2.0, 2.0)}
    return 8


def main() -> None:
    coordinate_checks = check_canonical_coordinates()
    gluing_checks = check_general_gluing_fibres()
    triangle_checks = check_triangle_classes()
    print(
        "orientation-carrier independent audit: PASS",
        f"coordinate_supports={coordinate_checks}",
        f"marginal_fibres={gluing_checks}",
        f"triangle_words={triangle_checks}",
    )


if __name__ == "__main__":
    main()
