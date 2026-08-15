#!/usr/bin/env python3
"""Finite audit of the rank-one moving-projection cut-code inequality.

The script evaluates the canonical k=0 moving-projection kernels from
Chapter 2 of OpenAI's ``Ten Advances in Mathematics and Theoretical Computer
Science`` on committed exact/good signing histograms.  It prints one JSON
object and writes no files.

Classification: floating-point finite evidence, not a certificate.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal


ROOT = Path(__file__).resolve().parents[1]


def load_json(relative: str) -> dict:
    with (ROOT / relative).open(encoding="utf-8") as handle:
        return json.load(handle)


def witness_matrix(n: int) -> list[list[int]]:
    if n <= 10:
        return load_json(f"computations/results/exact_m{n}.json")["matrix"]
    if n == 11:
        return load_json(
            "computations/results/nested_10_in_11_cap17.json"
        )["matrix"]
    if n == 12:
        return load_json(
            "computations/results/extension_nested_m11_to_12.json"
        )["parent_matrix"]
    if n == 13:
        matrix = load_json(
            "computations/results/conference_completion_m13.json"
        )["conference_matrix"]
        return [row[:-1] for row in matrix[:-1]]
    if n == 14:
        return load_json(
            "computations/results/heuristic_m14_from_conference.json"
        )["matrix"]
    raise ValueError(n)


def projective_energy_histogram(matrix: list[list[int]]) -> dict[int, int]:
    n = len(matrix)
    histogram: dict[int, int] = {}
    for mask in range(1 << (n - 1)):
        spin = [1] + [
            1 if (mask >> (i - 1)) & 1 else -1 for i in range(1, n)
        ]
        energy = sum(
            matrix[i][j] * spin[i] * spin[j]
            for i in range(n)
            for j in range(i + 1, n)
        )
        histogram[energy] = histogram.get(energy, 0) + 1
    return histogram


def rank_one_kernel_data(length: int, level: int) -> tuple[float, np.ndarray]:
    off_diagonal = np.array(
        [
            math.sqrt((i + 1) * (length - i)) / length
            for i in range(level)
        ]
    )
    eigenvalues, eigenvectors = eigh_tridiagonal(
        np.zeros(level + 1), off_diagonal
    )
    vector = eigenvectors[:, -1]
    if vector.sum() < 0:
        vector = -vector

    # In the paper w_i=sqrt(binomial(length,i))*v_i, and w_i/sum_j w_j
    # is the squared block amplitude of the rank-one moving vector.
    log_weight = np.array(
        [
            0.5
            * (
                math.lgamma(length + 1)
                - math.lgamma(i + 1)
                - math.lgamma(length - i + 1)
            )
            for i in range(level + 1)
        ]
    ) + np.log(np.maximum(vector, 1e-300))
    weight = np.exp(log_weight - log_weight.max())
    weight /= weight.sum()
    return float(eigenvalues[-1]), weight


def normalized_krawtchouk_sum(
    length: int, weight: np.ndarray, correlation: float
) -> float:
    total = float(weight[0])
    if len(weight) == 1:
        return total
    previous = 1.0
    current = correlation
    total += float(weight[1]) * current
    for i in range(1, len(weight) - 1):
        following = (
            length * correlation * current - i * previous
        ) / (length - i)
        total += float(weight[i + 1]) * following
        previous, current = current, following
    return total


def internal_remainder_average(
    n: int, length: int, weight: np.ndarray, eigenvalue: float
) -> float:
    total = 0.0
    for r in range(n + 1):
        q_r = ((n - 2 * r) ** 2 - n) / (n * (n - 1))
        for correlation in (q_r, -q_r):
            value = normalized_krawtchouk_sum(
                length, weight, correlation
            )
            total += (
                0.5
                * math.comb(n, r)
                * (correlation - eigenvalue)
                * value**2
            )
    return total / (2**n)


def coset_kernel_average(
    length: int, weight: np.ndarray, histogram: dict[int, int]
) -> float:
    projective_count = sum(histogram.values())
    total = 0.0
    for energy, count in histogram.items():
        correlation = energy / length
        positive = normalized_krawtchouk_sum(
            length, weight, correlation
        )
        negative = normalized_krawtchouk_sum(
            length, weight, -correlation
        )
        total += count * 0.5 * (positive**2 + negative**2)
    return total / projective_count


def audit_one(n: int, histogram: dict[int, int]) -> dict:
    length = math.comb(n, 2)
    cap = max(abs(energy) for energy in histogram)
    best: dict | None = None
    first_above_cap: dict | None = None
    final_level = min(length - 1, 3 * n)
    for level in range(1, final_level + 1):
        eigenvalue, weight = rank_one_kernel_data(length, level)
        remainder = internal_remainder_average(
            n, length, weight, eigenvalue
        )
        denominator = coset_kernel_average(length, weight, histogram)
        bound = eigenvalue - remainder / denominator
        candidate = {
            "level": level,
            "lambda": eigenvalue,
            "J_over_code_size": remainder,
            "T_over_code_size": denominator,
            "normalized_bound": bound,
            "energy_bound": length * bound,
        }
        if best is None or bound > best["normalized_bound"]:
            best = candidate
        if first_above_cap is None and eigenvalue > cap / length:
            first_above_cap = candidate
    assert best is not None
    return {
        "n": n,
        "edge_length": length,
        "actual_cap": cap,
        "actual_normalized_cap": cap / length,
        "levels_tested": [1, final_level],
        "best_rank_one_moving_bound": best,
        "first_level_with_lambda_above_actual_cap": first_above_cap,
        "rms_energy_bound": math.sqrt(length),
    }


def main() -> None:
    records = []
    for n in (6, 8, 10, 12, 14):
        records.append(
            audit_one(n, projective_energy_histogram(witness_matrix(n)))
        )

    conference_18 = load_json(
        "computations/results/conference_double_p17.json"
    )["conference_profile"]
    histogram_18 = {
        int(energy): int(count)
        for energy, count in conference_18["energy_histogram"].items()
    }
    records.append(audit_one(18, histogram_18))
    print(
        json.dumps(
            {"classification": "finite_numerical_audit", "records": records},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
