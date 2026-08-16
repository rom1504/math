#!/usr/bin/env python3
"""Exact finite checks for the metric-isometry bottleneck algebra.

The script checks, on several metrics of orders 2 through 5,

    D_(lambda,g) star D_(mu,h) = D_(min(lambda,mu), h o g),

the directed row-gap formula, full chain composition, the cumulative
entrywise perturbation bound, and the anisotropic projective-Hamming
coordinatewise-wedge law.  All arithmetic is integral and the random
families use the fixed seed printed in the output.
"""

from __future__ import annotations

import itertools
import json
import random
from collections.abc import Iterable, Sequence


Matrix = list[list[int]]
Permutation = tuple[int, ...]


def minplus(left: Matrix, right: Matrix) -> Matrix:
    size = len(left)
    return [
        [
            min(left[a][u] + right[u][t] for u in range(size))
            for t in range(size)
        ]
        for a in range(size)
    ]


def metric_kernel(
    metric: Matrix, isometry: Permutation, strength: int
) -> Matrix:
    size = len(metric)
    return [
        [strength * metric[t][isometry[a]] for t in range(size)]
        for a in range(size)
    ]


def compose_isometries(first: Permutation, second: Permutation) -> Permutation:
    """Return second o first, matching left-to-right kernel composition."""

    return tuple(second[first[a]] for a in range(len(first)))


def self_isometries(metric: Matrix) -> list[Permutation]:
    size = len(metric)
    return [
        permutation
        for permutation in itertools.permutations(range(size))
        if all(
            metric[permutation[a]][permutation[b]] == metric[a][b]
            for a in range(size)
            for b in range(size)
        )
    ]


def shortest_path_metric(size: int, rng: random.Random) -> Matrix:
    metric = [[0 for _ in range(size)] for _ in range(size)]
    for a in range(size):
        for b in range(a + 1, size):
            metric[a][b] = metric[b][a] = rng.randint(1, 9)
    for middle in range(size):
        for a in range(size):
            for b in range(size):
                metric[a][b] = min(
                    metric[a][b], metric[a][middle] + metric[middle][b]
                )
    return metric


def metric_families(size: int, rng: random.Random) -> list[Matrix]:
    candidates = [
        [[0 if a == b else 1 for b in range(size)] for a in range(size)],
        [[abs(a - b) for b in range(size)] for a in range(size)],
        [
            [min((a - b) % size, (b - a) % size) for b in range(size)]
            for a in range(size)
        ],
    ]
    candidates.extend(shortest_path_metric(size, rng) for _ in range(10))

    answer = []
    seen: set[tuple[tuple[int, ...], ...]] = set()
    for metric in candidates:
        key = tuple(tuple(row) for row in metric)
        if key not in seen:
            seen.add(key)
            answer.append(metric)
    return answer


def deterministic_sample(items: Sequence[Permutation], limit: int) -> list[Permutation]:
    if len(items) <= limit:
        return list(items)
    # Include both endpoints of the lexicographic list and representatives
    # spread deterministically throughout it.
    indices = {(i * (len(items) - 1)) // (limit - 1) for i in range(limit)}
    return [items[i] for i in sorted(indices)]


def max_entry_error(first: Matrix, second: Matrix) -> int:
    return max(
        abs(first[a][t] - second[a][t])
        for a in range(len(first))
        for t in range(len(first))
    )


def directed_gap(matrix: Matrix, first_row: int, second_row: int) -> int:
    return max(
        matrix[first_row][t] - matrix[second_row][t]
        for t in range(len(matrix))
    )


def verify_exact_two_factor(metrics: Iterable[Matrix]) -> tuple[int, int]:
    composition_checks = 0
    gap_checks = 0
    for metric in metrics:
        size = len(metric)
        isometries = self_isometries(metric)
        sample = deterministic_sample(isometries, 24)
        for first in sample:
            for second in sample:
                composed = compose_isometries(first, second)
                for first_strength in (0, 1, 3, 7):
                    for second_strength in (0, 2, 5):
                        observed = minplus(
                            metric_kernel(metric, first, first_strength),
                            metric_kernel(metric, second, second_strength),
                        )
                        predicted = metric_kernel(
                            metric,
                            composed,
                            min(first_strength, second_strength),
                        )
                        assert observed == predicted
                        composition_checks += 1

        for isometry in sample:
            for strength in (0, 1, 4):
                kernel = metric_kernel(metric, isometry, strength)
                for first_row in range(size):
                    for second_row in range(size):
                        assert directed_gap(kernel, first_row, second_row) == (
                            strength * metric[first_row][second_row]
                        )
                        gap_checks += 1
    return composition_checks, gap_checks


def perturb_matrix(
    matrix: Matrix, radius: int, rng: random.Random
) -> Matrix:
    return [
        [entry + rng.randint(-radius, radius) for entry in row]
        for row in matrix
    ]


def verify_chains(
    metrics: Sequence[Matrix], rng: random.Random
) -> tuple[int, int, int]:
    exact_chain_checks = 0
    uniform_error_checks = 0
    directed_error_checks = 0

    for metric in metrics:
        size = len(metric)
        isometries = self_isometries(metric)
        for length in range(1, 9):
            for _ in range(30):
                labels = [rng.choice(isometries) for _ in range(length)]
                strengths = [rng.randint(0, 12) for _ in range(length)]
                radii = [rng.randint(0, 3) for _ in range(length)]

                exact = metric_kernel(metric, labels[0], strengths[0])
                noisy = perturb_matrix(exact, radii[0], rng)
                holonomy = labels[0]
                for label, strength, radius in zip(
                    labels[1:], strengths[1:], radii[1:]
                ):
                    factor = metric_kernel(metric, label, strength)
                    exact = minplus(exact, factor)
                    noisy = minplus(noisy, perturb_matrix(factor, radius, rng))
                    holonomy = compose_isometries(holonomy, label)

                bottleneck = min(strengths)
                predicted = metric_kernel(metric, holonomy, bottleneck)
                assert exact == predicted
                exact_chain_checks += 1

                total_radius = sum(radii)
                assert max_entry_error(noisy, predicted) <= total_radius
                uniform_error_checks += 1

                for first_row in range(size):
                    for second_row in range(size):
                        ideal_gap = bottleneck * metric[first_row][second_row]
                        assert abs(
                            directed_gap(noisy, first_row, second_row) - ideal_gap
                        ) <= 2 * total_radius
                        directed_error_checks += 1

    return exact_chain_checks, uniform_error_checks, directed_error_checks


def canonical_projective_word(word: tuple[int, ...]) -> tuple[int, ...]:
    complement = tuple(1 - bit for bit in word)
    return min(word, complement)


def projective_words(width: int) -> list[tuple[int, ...]]:
    return sorted(
        {
            canonical_projective_word(word)
            for word in itertools.product((0, 1), repeat=width)
        }
    )


def projective_hamming_distance(
    first: tuple[int, ...], second: tuple[int, ...], weights: Sequence[int]
) -> int:
    oriented = sum(
        weight
        for bit_first, bit_second, weight in zip(first, second, weights)
        if bit_first != bit_second
    )
    return min(oriented, sum(weights) - oriented)


def act_on_coordinates(
    permutation: Permutation, word: tuple[int, ...]
) -> tuple[int, ...]:
    # The coordinate originally at i moves to permutation[i].
    image = [0 for _ in word]
    for old_coordinate, new_coordinate in enumerate(permutation):
        image[new_coordinate] = word[old_coordinate]
    return canonical_projective_word(tuple(image))


def pushforward_weights(
    permutation: Permutation, weights: Sequence[int]
) -> list[int]:
    image = [0 for _ in weights]
    for old_coordinate, new_coordinate in enumerate(permutation):
        image[new_coordinate] = weights[old_coordinate]
    return image


def verify_anisotropic_projective_hamming(
    rng: random.Random,
) -> int:
    """Check B_(ell,g) B_(m,h) = B_(m wedge h_*ell,hg)."""

    checks = 0
    for width in range(2, 7):
        words = projective_words(width)
        for _ in range(30):
            first_permutation = tuple(rng.sample(range(width), width))
            second_permutation = tuple(rng.sample(range(width), width))
            first_weights = [rng.randint(0, 9) for _ in range(width)]
            second_weights = [rng.randint(0, 9) for _ in range(width)]

            transported = pushforward_weights(
                second_permutation, first_weights
            )
            bottleneck_weights = [
                min(first, second)
                for first, second in zip(transported, second_weights)
            ]
            holonomy = compose_isometries(
                first_permutation, second_permutation
            )

            for first_word in words:
                first_center = act_on_coordinates(
                    first_permutation, first_word
                )
                predicted_center = act_on_coordinates(holonomy, first_word)
                for last_word in words:
                    observed = min(
                        projective_hamming_distance(
                            intermediate, first_center, first_weights
                        )
                        + projective_hamming_distance(
                            last_word,
                            act_on_coordinates(second_permutation, intermediate),
                            second_weights,
                        )
                        for intermediate in words
                    )
                    predicted = projective_hamming_distance(
                        last_word, predicted_center, bottleneck_weights
                    )
                    assert observed == predicted
                    checks += 1
    return checks


def main() -> None:
    seed = 20260816
    rng = random.Random(seed)
    metrics = [
        metric
        for size in range(2, 6)
        for metric in metric_families(size, rng)
    ]

    composition_checks, gap_checks = verify_exact_two_factor(metrics)
    exact_chains, uniform_errors, directed_errors = verify_chains(metrics, rng)
    anisotropic_checks = verify_anisotropic_projective_hamming(rng)

    print(
        json.dumps(
            {
                "seed": seed,
                "metric_instances": len(metrics),
                "two_factor_composition_checks": composition_checks,
                "directed_gap_checks": gap_checks,
                "exact_chain_checks": exact_chains,
                "uniform_perturbation_checks": uniform_errors,
                "perturbed_directed_gap_checks": directed_errors,
                "anisotropic_projective_hamming_checks": anisotropic_checks,
            },
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
