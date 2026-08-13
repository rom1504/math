#!/usr/bin/env python3
"""Independently verify the global Fourier-layer formulas by exact arithmetic."""

from __future__ import annotations

import itertools
import json
import math
import random

import numpy as np


def choose(n: int, k: int) -> int:
    return math.comb(n, k) if 0 <= k <= n else 0


def operator(
    matrix: np.ndarray, levels: tuple[int, ...]
) -> tuple[np.ndarray, list[tuple[int, ...]]]:
    n = len(matrix)
    subsets = [
        subset
        for level in levels
        for subset in itertools.combinations(range(n), level)
    ]
    position = {subset: index for index, subset in enumerate(subsets)}
    result = np.zeros((len(subsets), len(subsets)), dtype=np.int64)
    for row, subset_tuple in enumerate(subsets):
        subset = set(subset_tuple)
        for i in range(n):
            for j in range(i + 1, n):
                target = tuple(sorted(subset.symmetric_difference({i, j})))
                if target in position:
                    result[row, position[target]] = int(matrix[i, j])
    return result, subsets


def energy(matrix: np.ndarray, spin: np.ndarray) -> int:
    return int(spin @ matrix @ spin) // 2


def cycle_sum(matrix: np.ndarray) -> int:
    result = 0
    for quadruple in itertools.combinations(range(len(matrix)), 4):
        i, j, k, ell = quadruple
        result += int(matrix[i, j] * matrix[j, k] * matrix[k, ell] * matrix[ell, i])
        result += int(matrix[i, j] * matrix[j, ell] * matrix[ell, k] * matrix[k, i])
        result += int(matrix[i, k] * matrix[k, j] * matrix[j, ell] * matrix[ell, i])
    return result


def trace_four_formula(matrix: np.ndarray, level: int) -> int:
    n = len(matrix)
    coefficient = (
        32 * choose(n - 4, level - 2)
        + 8 * choose(n - 4, level - 1)
        + 8 * choose(n - 4, level - 3)
    )
    johnson_trace = 0
    for j in range(min(level, n - level) + 1):
        multiplicity = choose(n, j) - choose(n, j - 1)
        eigenvalue = (level - j) * (n - level - j) - j
        johnson_trace += multiplicity * eigenvalue**4
    return (
        johnson_trace
        - 3 * coefficient * choose(n, 4)
        + coefficient * cycle_sum(matrix)
    )


def band_rho_numerator(n: int, levels: tuple[int, ...]) -> tuple[int, int]:
    dimension = sum(choose(n, level) for level in levels)
    boundary = choose(n - 2, levels[0] - 2) + choose(n - 2, levels[-1])
    return dimension - boundary, dimension


def main() -> None:
    rng = random.Random(20260813)
    trace_cases = 0
    band_cases = 0
    for n in range(4, 10):
        for _sample in range(4):
            matrix = np.zeros((n, n), dtype=np.int64)
            for i in range(n):
                for j in range(i + 1, n):
                    matrix[i, j] = matrix[j, i] = rng.choice((-1, 1))

            for level in range(n + 1):
                layer, _subsets = operator(matrix, (level,))
                actual = int(np.trace(layer @ layer @ layer @ layer))
                predicted = trace_four_formula(matrix, level)
                if actual != predicted:
                    raise AssertionError((n, level, actual, predicted))
                trace_cases += 1

            for parity in (0, 1):
                available = [level for level in range(n + 1) if level % 2 == parity]
                for low in range(len(available)):
                    for high in range(low, len(available)):
                        levels = tuple(available[low : high + 1])
                        layer, subsets = operator(matrix, levels)
                        numerator, _dimension = band_rho_numerator(n, levels)
                        for _spin_sample in range(3):
                            spin = np.asarray(
                                [rng.choice((-1, 1)) for _ in range(n)],
                                dtype=np.int64,
                            )
                            vector = np.asarray(
                                [math.prod(int(spin[i]) for i in subset) for subset in subsets],
                                dtype=np.int64,
                            )
                            actual = int(vector @ layer @ vector)
                            predicted = numerator * energy(matrix, spin)
                            if actual != predicted:
                                raise AssertionError(
                                    (n, levels, actual, predicted)
                                )
                        band_cases += 1

    output = {
        "band_rho_cases": band_cases,
        "orders": "4 through 9",
        "schema": "global-layer-formulas-independent-check-v1",
        "seed": 20260813,
        "spin_checks_per_band_case": 3,
        "status": "all exact integer identities passed",
        "trace_four_cases": trace_cases,
    }
    print(json.dumps(output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
