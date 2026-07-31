#!/usr/bin/env python3
"""Compute the exact orbit-Fourier transform behind the phi_6 profile."""

from __future__ import annotations

import argparse
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np

from phase2_restriction_state_audit import (
    class_map,
    root_gauge_code,
    signing_from_code,
)


def popcount(value: int) -> int:
    return bin(value).count("1")


def invariant_edge_orbits(
    n: int,
) -> tuple[list[tuple[tuple[int, int], ...]], list[tuple[int, ...]]]:
    edges = list(itertools.combinations(range(n), 2))
    edge_index = {edge: index for index, edge in enumerate(edges)}
    permutations = list(itertools.permutations(range(n)))

    def move(mask: int, permutation: tuple[int, ...]) -> int:
        result = 0
        for index, (i, j) in enumerate(edges):
            if mask & (1 << index):
                image = tuple(sorted((permutation[i], permutation[j])))
                result |= 1 << edge_index[image]
        return result

    by_canonical: dict[int, tuple[int, ...]] = {}
    for mask in range(1 << len(edges)):
        if popcount(mask) % 2:
            continue
        degrees = [0] * n
        for index, (i, j) in enumerate(edges):
            if mask & (1 << index):
                degrees[i] += 1
                degrees[j] += 1
        if any(degree % 2 for degree in degrees):
            continue
        orbit = tuple(sorted({move(mask, permutation) for permutation in permutations}))
        by_canonical[orbit[0]] = orbit
    ordered = [
        by_canonical[key]
        for key in sorted(by_canonical, key=lambda item: (popcount(item), item))
    ]
    return edges, ordered


def graph_descriptor(
    n: int, edges: list[tuple[int, int]], mask: int, orbit_size: int
) -> dict[str, object]:
    chosen = [edge for index, edge in enumerate(edges) if mask & (1 << index)]
    degrees = [0] * n
    adjacency = [set() for _ in range(n)]
    for i, j in chosen:
        degrees[i] += 1
        degrees[j] += 1
        adjacency[i].add(j)
        adjacency[j].add(i)
    seen: set[int] = set()
    component_sizes = []
    for start in range(n):
        if start in seen:
            continue
        seen.add(start)
        stack = [start]
        component = []
        while stack:
            vertex = stack.pop()
            component.append(vertex)
            for neighbor in adjacency[vertex]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
        component_sizes.append(len(component))
    triangles = sum(
        1
        for i, j, k in itertools.combinations(range(n), 3)
        if j in adjacency[i] and k in adjacency[i] and k in adjacency[j]
    )
    return {
        "canonical_edges": [list(edge) for edge in chosen],
        "edge_count": len(chosen),
        "degree_sequence": sorted(degrees, reverse=True),
        "component_sizes": sorted(component_sizes, reverse=True),
        "triangle_count": triangles,
        "labeled_orbit_size": orbit_size,
    }


def character_transform(n: int) -> dict[str, object]:
    labels, class_count = class_map(n)
    representative_codes = [labels.index(label) for label in range(class_count)]
    edges, graph_orbits = invariant_edge_orbits(n)
    matrix = []
    for orbit in graph_orbits:
        row = []
        for code in representative_codes:
            signing = signing_from_code(code, n)
            value = 0
            for mask in orbit:
                product = 1
                for index, (i, j) in enumerate(edges):
                    if mask & (1 << index):
                        product *= int(signing[i, j])
                value += product
            row.append(value)
        matrix.append(row)
    if len(matrix) != class_count:
        raise AssertionError((n, len(matrix), class_count))
    work = [[Fraction(value) for value in row] for row in matrix]
    rank = 0
    determinant = Fraction(1)
    for column in range(class_count):
        pivot = next(
            (row for row in range(rank, class_count) if work[row][column]), None
        )
        if pivot is None:
            continue
        if pivot != rank:
            work[rank], work[pivot] = work[pivot], work[rank]
            determinant *= -1
        pivot_value = work[rank][column]
        determinant *= pivot_value
        work[rank] = [value / pivot_value for value in work[rank]]
        for row in range(class_count):
            if row == rank or not work[row][column]:
                continue
            factor = work[row][column]
            work[row] = [
                left - factor * right for left, right in zip(work[row], work[rank])
            ]
        rank += 1
    return {
        "order": n,
        "switching_permutation_global_negation_classes": class_count,
        "representative_root_codes": representative_codes,
        "invariant_graph_orbits": [
            graph_descriptor(n, edges, orbit[0], len(orbit)) for orbit in graph_orbits
        ],
        "character_transform_rows_graphs_columns_signing_classes": matrix,
        "rank_over_Q": rank,
        "determinant": str(determinant),
    }


def six_deck_incidence() -> dict[str, object]:
    labels6, class_count6 = class_map(6)
    representatives6 = [labels6.index(label) for label in range(class_count6)]
    result: dict[str, object] = {}
    for size in (4, 5):
        labels, class_count = class_map(size)
        incidence = [[0] * class_count6 for _ in range(class_count)]
        for column, code in enumerate(representatives6):
            matrix = signing_from_code(code, 6)
            for vertices in itertools.combinations(range(6), size):
                child = matrix[np.ix_(vertices, vertices)]
                row = labels[root_gauge_code(child)]
                incidence[row][column] += 1
        if any(sum(incidence[row][column] for row in range(class_count))
               != math.comb(6, size) for column in range(class_count6)):
            raise AssertionError((size, incidence))
        result[str(size)] = {
            "rows": f"order-{size} classes",
            "columns": "order-6 classes",
            "incidence_matrix": incidence,
            "recovery_formula": (
                f"h_{size} = incidence*h_6/binom(n-{size},{6-size})"
            ),
        }
    return result


def orbit_moment(
    matrix: np.ndarray, size: int, canonical_edges: list[list[int]]
) -> int:
    local_edges = list(itertools.combinations(range(size), 2))
    local_index = {edge: index for index, edge in enumerate(local_edges)}
    base_mask = 0
    for edge_list in canonical_edges:
        edge = tuple(edge_list)
        base_mask |= 1 << local_index[edge]
    permutations = itertools.permutations(range(size))
    masks = set()
    for permutation in permutations:
        moved = 0
        for index, (i, j) in enumerate(local_edges):
            if base_mask & (1 << index):
                image = tuple(sorted((permutation[i], permutation[j])))
                moved |= 1 << local_index[image]
        masks.add(moved)
    result = 0
    for vertices in itertools.combinations(range(len(matrix)), size):
        for mask in masks:
            product = 1
            for index, (i, j) in enumerate(local_edges):
                if mask & (1 << index):
                    product *= int(matrix[vertices[i], vertices[j]])
            result += product
    return result


def verify_moment_formulas(transforms: dict[int, dict[str, object]]) -> list[dict[str, object]]:
    records = []
    for n in (6, 7, 8, 9):
        generator = np.random.default_rng(20260801 + n)
        matrix = np.zeros((n, n), dtype=np.int64)
        for i in range(n):
            for j in range(i + 1, n):
                matrix[i, j] = matrix[j, i] = generator.choice((-1, 1))
        moments: dict[tuple[int, int], int] = {}
        for size, transform in transforms.items():
            for index, descriptor in enumerate(transform["invariant_graph_orbits"]):
                moments[(size, index)] = orbit_moment(
                    matrix, size, descriptor["canonical_edges"]
                )
        edge_count = math.comb(n, 2)
        z4 = moments[(4, 1)]
        trace4_formula = n * (n - 1) * (2 * n - 3) + 8 * z4
        trace6_formula = (
            2 * math.comb(n, 2)
            + 60 * math.comb(n, 3)
            + 120 * math.comb(n, 4)
            + 120 * moments[(4, 1)]
            + 48 * moments[(5, 1)]
            + 24 * moments[(5, 2)]
            + 12 * moments[(6, 3)]
        )
        energy4_formula = 3 * edge_count**2 - 2 * edge_count + 24 * z4
        energy6_formula = (
            edge_count
            + 15 * edge_count * (edge_count - 1)
            + 90 * math.comb(edge_count, 3)
            + (360 * edge_count - 960) * z4
            + 720
            * (moments[(5, 2)] + moments[(6, 3)] + moments[(6, 4)])
        )
        trace4 = int(np.trace(np.linalg.matrix_power(matrix, 4)))
        trace6 = int(np.trace(np.linalg.matrix_power(matrix, 6)))
        spins = np.asarray(
            [
                [1 if not (code >> i) & 1 else -1 for i in range(n)]
                for code in range(1 << n)
            ],
            dtype=np.int64,
        )
        energies = np.einsum(
            "bi,ij,bj->b", spins, matrix, spins, optimize=True
        ) // 2
        energy4 = sum(int(value) ** 4 for value in energies) // len(energies)
        energy6 = sum(int(value) ** 6 for value in energies) // len(energies)
        checks = {
            "trace4": trace4 == trace4_formula,
            "trace6": trace6 == trace6_formula,
            "uniform_energy4": energy4 == energy4_formula,
            "uniform_energy6": energy6 == energy6_formula,
        }
        if not all(checks.values()):
            raise AssertionError((n, checks))
        records.append(
            {
                "order": n,
                "matrix_seed": 20260801 + n,
                "checks": checks,
                "trace4": trace4,
                "trace6": trace6,
                "uniform_energy_fourth_moment": energy4,
                "uniform_energy_sixth_moment": energy6,
            }
        )
    return records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    transforms = {size: character_transform(size) for size in (4, 5, 6)}
    output = {
        "schema": "quadratic-signing-phi6-moment-transform-v1",
        "classification": "exact finite invariant theory and arithmetic verification",
        "normalization": (
            "orbit moment sums, over every vertex subset, the permutation-orbit "
            "of an even-edge Eulerian sign monomial"
        ),
        "transforms": {str(size): value for size, value in transforms.items()},
        "six_deck_determines_lower_profiles": six_deck_incidence(),
        "formula_verifications": verify_moment_formulas(transforms),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps(output, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
