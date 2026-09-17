#!/usr/bin/env python3
"""Search general linear-code deep holes for nonmonotone noise curves.

This can falsify an over-general convexity/code-distance lemma, but says nothing
by itself about complete-graph cut codes. Every code here contains all ones,
so its support-function norm is invariant under global coefficient reversal.
"""

from __future__ import annotations

import argparse
from collections import deque
import json
from pathlib import Path

import numpy as np
import sympy as sp


def fwht(values: np.ndarray) -> np.ndarray:
    result = values.copy()
    step = 1
    while step < len(result):
        blocks = result.reshape(-1, 2 * step)
        a = blocks[:, :step].copy()
        b = blocks[:, step:].copy()
        blocks[:, :step] = a + b
        blocks[:, step:] = a - b
        step *= 2
    return result


def examine(columns: np.ndarray, rank: int, min_dual_weight: int) -> dict | None:
    dimension = len(columns)
    syndromes = np.arange(1 << rank, dtype=np.int64)
    parity = np.array([bin(i).count("1") % 2 for i in syndromes], dtype=np.int64)
    weights = np.zeros(len(syndromes), dtype=np.int64)
    for col in columns:
        weights += parity[syndromes & col]
    if np.any(weights % 2):
        raise AssertionError("global all-ones codeword missing")
    if weights[1:].min() < min_dual_weight:
        return None
    distance = np.full(len(syndromes), dimension + 1, dtype=np.int64)
    distance[0] = 0
    queue = deque([0])
    while queue:
        vertex = queue.popleft()
        for col in columns:
            neighbor = vertex ^ int(col)
            if distance[neighbor] > distance[vertex] + 1:
                distance[neighbor] = distance[vertex] + 1
                queue.append(neighbor)
    if distance.max() > dimension:
        return None
    value = dimension - 2 * distance
    transform = fwht(value)
    levels = np.zeros((dimension + 1, len(syndromes)), dtype=np.int64)
    for k in range(0, dimension + 1, 2):
        levels[k] = fwht(np.where(weights == k, transform, 0))
    deepest = np.flatnonzero(distance == distance.max())
    grid = np.concatenate((np.linspace(0.001, 0.1, 25), np.linspace(0.1, 0.999, 150)))
    coefficients = levels[:, deepest]
    for t in grid:
        derivative = np.sum(np.arange(dimension + 1)[:, None] * coefficients *
                            np.power(t, np.maximum(np.arange(dimension + 1) - 1, 0))[:, None], axis=0)
        candidates = np.flatnonzero(derivative > 1e-8)
        if len(candidates):
            idx = int(candidates[np.argmax(derivative[candidates])])
            syndrome = int(deepest[idx])
            numerators = coefficients[:, idx]
            rational_t = sp.Rational(str(float(t)))
            exact_derivative = sum(int(numerators[k]) * k * rational_t ** (k - 1)
                                   for k in range(1, dimension + 1)) / len(syndromes)
            assert exact_derivative > 0
            return {"ambient_dimension": dimension, "parity_check_rank": rank,
                    "parity_check_columns_integer": columns.tolist(),
                    "deepest_syndrome": syndrome, "covering_radius": int(distance.max()),
                    "minimum_norm": int(value[syndrome]),
                    "minimum_dual_weight": int(weights[1:].min()),
                    "noise_polynomial_numerators": numerators.tolist(),
                    "noise_polynomial_denominator": len(syndromes),
                    "positive_derivative_at_t": str(rational_t),
                    "positive_derivative_exact": str(exact_derivative),
                    "positive_derivative_float": float(exact_derivative)}
    return None


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seed", type=int, default=20260917)
    parser.add_argument("--trials", type=int, default=10000)
    parser.add_argument("--min-dual-weight", type=int, default=4)
    parser.add_argument("--output", type=Path,
                        default=Path("tmp/paper_portfolio_2026_09_17/localization/code_noise_search.json"))
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    found = None
    for trial in range(args.trials):
        rank = int(rng.integers(3, 9))
        dimension = int(rng.integers(rank + 2, 3 * rank + 1))
        columns = rng.integers(1, 1 << rank, dimension - 1, dtype=np.int64)
        columns = np.concatenate((columns, [np.bitwise_xor.reduce(columns)]))
        found = examine(columns, rank, args.min_dual_weight)
        if found:
            break
    result = {"seed": args.seed, "trials_completed": trial + 1, "counterexample": found}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
