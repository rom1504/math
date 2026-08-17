#!/usr/bin/env python3
"""Finite verifier for the multiscale partition affine-shell theorem.

The default run exhausts all hollow signings through order five and checks
random signings at larger orders.  For every q it verifies the two one-sided
partition budgets, the subset-completion bound, the oriented affine cube,
and both the absolute and oriented one-sided star-frame responses.
"""

from __future__ import annotations

import argparse
import itertools
import random


def spins(n: int):
    return itertools.product((-1, 1), repeat=n)


def energy(matrix, spin, vertices=None) -> int:
    if vertices is None:
        vertices = range(len(matrix))
    vertices = list(vertices)
    return sum(
        matrix[i][j] * spin[i] * spin[j]
        for position, i in enumerate(vertices)
        for j in vertices[position + 1 :]
    )


def restricted_caps(matrix, vertices) -> tuple[int, int]:
    vertices = list(vertices)
    values = []
    for local_spin in spins(len(vertices)):
        spin = [1] * len(matrix)
        for i, value in zip(vertices, local_spin):
            spin[i] = value
        values.append(energy(matrix, spin, vertices))
    return max(values), -min(values)


def flip(spin, vertices):
    result = list(spin)
    for i in vertices:
        result[i] *= -1
    return tuple(result)


def verify(matrix) -> None:
    n = len(matrix)
    landscape = [(tuple(x), energy(matrix, x)) for x in spins(n)]
    ground, ground_energy = max(landscape, key=lambda item: abs(item[1]))
    cap = abs(ground_energy)
    orientation = 1 if ground_energy >= 0 else -1

    switched = [
        [orientation * ground[i] * matrix[i][j] * ground[j] for j in range(n)]
        for i in range(n)
    ]
    local_fields = [sum(row) for row in switched]
    assert all(value >= 0 for value in local_fields)
    assert sum(local_fields) == 2 * cap

    for block_count in range(2, n + 1):
        blocks = [
            list(range(a * n // block_count, (a + 1) * n // block_count))
            for a in range(block_count)
        ]
        block_data = []
        positive_sum = 0
        negative_sum = 0
        for block in blocks:
            positive, negative = restricted_caps(switched, block)
            positive_sum += positive
            negative_sum += negative
            field_mass = sum(local_fields[i] for i in block)
            block_data.append((2 * field_mass + 4 * negative, block, negative))

            for mask in range(1 << len(block)):
                subset = [block[j] for j in range(len(block)) if mask >> j & 1]
                subset_mass = sum(
                    switched[i][j]
                    for position, i in enumerate(subset)
                    for j in subset[position + 1 :]
                )
                assert -subset_mass <= negative

        assert positive_sum <= cap
        assert negative_sum <= cap

        score, block, _ = min(block_data)
        assert score * block_count <= 8 * cap
        assert len(block) >= n // block_count

        cube = []
        for mask in range(1 << len(block)):
            subset = [block[j] for j in range(len(block)) if mask >> j & 1]
            member = flip(ground, subset)
            cube.append(member)
            assert orientation * energy(matrix, member) >= cap - score
            assert orientation * energy(matrix, member) >= cap - 8 * cap / block_count

        representatives = {min(x, tuple(-v for v in x)) for x in cube}
        assert len(representatives) == 1 << len(block)

        # Largest even subframe, with its abstract endpoint language.
        even_size = len(block) - len(block) % 2
        interface = block[:even_size]
        ports = [ground] + [flip(ground, [i]) for i in interface]
        for endpoint in spins(len(ports)):
            field = [
                sum(sign * port[i] for sign, port in zip(endpoint, ports))
                for i in range(n)
            ]
            assert all(value != 0 for value in field)
            selector = tuple(1 if value > 0 else -1 for value in field)
            selector_orbit = {
                member for member in cube
            } | {
                tuple(-v for v in member) for member in cube
            }
            assert selector in selector_orbit

            field_norm = sum(abs(value) for value in field)
            absolute_response = max(
                abs(value) + sum(g * y[i] for i, g in enumerate(field))
                for y, value in landscape
            )
            oriented_response = max(
                orientation * value + sum(g * y[i] for i, g in enumerate(field))
                for y, value in landscape
            )
            for response in (absolute_response, oriented_response):
                gap = cap + field_norm - response
                assert -1e-12 <= gap <= 8 * cap / block_count + 1e-12


def signing_matrix(n: int, edge_signs):
    matrix = [[0] * n for _ in range(n)]
    cursor = 0
    for i in range(n):
        for j in range(i + 1, n):
            matrix[i][j] = matrix[j][i] = edge_signs[cursor]
            cursor += 1
    return matrix


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=20260817)
    parser.add_argument("--random-per-order", type=int, default=12)
    parser.add_argument("--max-random-order", type=int, default=9)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    checked = 0

    for n in range(2, 6):
        edge_count = n * (n - 1) // 2
        for edge_signs in spins(edge_count):
            verify(signing_matrix(n, edge_signs))
            checked += 1

    for n in range(6, args.max_random_order + 1):
        edge_count = n * (n - 1) // 2
        for _ in range(args.random_per_order):
            edge_signs = tuple(rng.choice((-1, 1)) for _ in range(edge_count))
            verify(signing_matrix(n, edge_signs))
            checked += 1

    print(f"PASS: {checked} signings through order {args.max_random_order}")


if __name__ == "__main__":
    main()
