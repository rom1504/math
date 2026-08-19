#!/usr/bin/env python3
"""Exhaust the universal-double construction on small pressure minimizers.

For a hollow signing ``A`` of order ``r`` this audits

    K_0(A) = [[A, A],     [A, -A]]
    K_+(A) = [[A, A + I], [A + I, -A]].

``K_0`` has ``r`` zero matching edges and is only the fractional core;
``K_+`` is the required integral signing.  At scaled temperature ``beta``
the child and parent raw temperatures are respectively
``s=beta/sqrt(r)`` and ``t=beta/sqrt(2r)``.

Every switching orbit is represented by fixing the edges incident with
vertex zero to ``+1``.  The enumeration is exhaustive through ``--max-r``.
Only matrices attaining the minimum normalized cosh pressure at a requested
temperature are subsequently quotiented by permutation.  Global negation is
*not* quotiented out: child pressure does not see it, while ``K_+`` can.

Floating-point evaluations choose among a finite list of exact energy
histograms.  The output records the gap to the next histogram at every grid
point, as well as exact integer energy histograms for independent checking.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


DEFAULT_BETAS = [0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0, 16.0]


def projective_spins(n: int) -> np.ndarray:
    """All Boolean spins with the first coordinate fixed to +1."""
    if n == 0:
        return np.empty((1, 0), dtype=np.int8)
    result = np.ones((1 << max(n - 1, 0), n), dtype=np.int8)
    codes = np.arange(len(result), dtype=np.uint64)
    for index in range(1, n):
        result[:, index] = 1 - 2 * (
            (codes >> (index - 1)) & 1
        ).astype(np.int8)
    return result


def free_edges(n: int) -> list[tuple[int, int]]:
    return [(i, j) for i in range(1, n) for j in range(i + 1, n)]


def matrix_from_root_code(n: int, code: int) -> np.ndarray:
    matrix = np.ones((n, n), dtype=np.int8)
    np.fill_diagonal(matrix, 0)
    for bit, (i, j) in enumerate(free_edges(n)):
        if (code >> bit) & 1:
            matrix[i, j] = matrix[j, i] = -1
    return matrix


def energy_values(matrix: np.ndarray) -> np.ndarray:
    spins = projective_spins(len(matrix))
    return (
        np.einsum(
            "bi,ij,bj->b",
            spins,
            matrix.astype(np.int64),
            spins,
            optimize=True,
        )
        // 2
    )


def signed_histogram(values: np.ndarray) -> dict[str, int]:
    unique, counts = np.unique(values, return_counts=True)
    return {str(int(value)): int(count) for value, count in zip(unique, counts)}


def absolute_histogram_key(values: np.ndarray, edge_count: int) -> tuple[int, ...]:
    return tuple(
        int(value)
        for value in np.bincount(np.abs(values), minlength=edge_count + 1)
    )


def histogram_record(key: tuple[int, ...]) -> dict[str, int]:
    return {str(energy): count for energy, count in enumerate(key) if count}


def log_mean_cosh(values: np.ndarray, raw_temperature: float) -> float:
    scaled = raw_temperature * np.asarray(values, dtype=np.float64)
    peak = float(np.max(np.abs(scaled))) if len(scaled) else 0.0
    terms = 0.5 * (np.exp(scaled - peak) + np.exp(-scaled - peak))
    return peak + math.log(float(np.mean(terms)))


def pressure_from_absolute_histogram(
    key: tuple[int, ...], raw_temperature: float
) -> float:
    peak = raw_temperature * (len(key) - 1)
    total = math.fsum(
        count
        * 0.5
        * (
            math.exp(raw_temperature * energy - peak)
            + math.exp(-raw_temperature * energy - peak)
        )
        for energy, count in enumerate(key)
    )
    return peak + math.log(total / sum(key))


def enumerate_root_energy_data(
    n: int,
) -> tuple[list[np.ndarray], dict[tuple[int, ...], list[int]]]:
    """Return projective energies and histogram-to-root-code fibres."""
    edges = list(itertools.combinations(range(n), 2))
    variable_edges = free_edges(n)
    signing_count = 1 << len(variable_edges)
    spins = projective_spins(n)
    products = np.stack(
        [spins[:, i] * spins[:, j] for i, j in edges], axis=1
    ) if edges else np.empty((len(spins), 0), dtype=np.int8)

    codes = np.arange(signing_count, dtype=np.uint64)
    signs = np.ones((signing_count, len(edges)), dtype=np.int8)
    edge_position = {edge: index for index, edge in enumerate(edges)}
    for bit, edge in enumerate(variable_edges):
        signs[:, edge_position[edge]] = 1 - 2 * (
            (codes >> bit) & 1
        ).astype(np.int8)
    all_values = signs.astype(np.int16) @ products.T.astype(np.int16)

    values_by_code: list[np.ndarray] = []
    fibres: dict[tuple[int, ...], list[int]] = defaultdict(list)
    for code, values in enumerate(all_values):
        copied = values.astype(np.int64, copy=True)
        values_by_code.append(copied)
        fibres[absolute_histogram_key(copied, len(edges))].append(code)
    return values_by_code, dict(fibres)


def transformed_root_code(n: int, code: int, permutation: tuple[int, ...]) -> int:
    """Permute a root-gauged signing and switch it back to root gauge."""
    positions = {edge: bit for bit, edge in enumerate(free_edges(n))}

    def negative_bit(i: int, j: int) -> int:
        if i > j:
            i, j = j, i
        if i == 0:
            return 0
        return (code >> positions[(i, j)]) & 1

    result = 0
    old_root = permutation[0]
    for bit, (i, j) in enumerate(free_edges(n)):
        old_i = permutation[i]
        old_j = permutation[j]
        value = (
            negative_bit(old_i, old_j)
            ^ negative_bit(old_root, old_i)
            ^ negative_bit(old_root, old_j)
        )
        result |= value << bit
    return result


def permutation_classes(n: int, selected_codes: set[int]) -> list[dict[str, object]]:
    """Quotient selected root-gauge codes by vertex permutation only."""
    permutations = list(itertools.permutations(range(n)))
    remaining = set(selected_codes)
    records: list[dict[str, object]] = []
    while remaining:
        seed = min(remaining)
        orbit = {
            transformed_root_code(n, seed, permutation)
            for permutation in permutations
        }
        selected_members = sorted(remaining.intersection(orbit))
        remaining.difference_update(orbit)
        records.append(
            {
                "canonical_root_code": min(orbit),
                "selected_root_code_count": len(selected_members),
                "selected_root_code_examples": selected_members[:8],
                "full_root_gauge_permutation_orbit_size": len(orbit),
                "_selected_root_codes": selected_members,
            }
        )
    records.sort(key=lambda record: int(record["canonical_root_code"]))
    return records


def universal_double(matrix: np.ndarray, matching: int) -> np.ndarray:
    n = len(matrix)
    bridge = matrix.copy()
    if matching:
        bridge = bridge + matching * np.eye(n, dtype=np.int8)
    return np.block([[matrix, bridge], [bridge, -matrix]]).astype(np.int8)


def induced_flip(matrix: np.ndarray, subset_mask: int) -> np.ndarray:
    result = matrix.copy()
    n = len(matrix)
    for i in range(n):
        if not ((subset_mask >> i) & 1):
            continue
        for j in range(i + 1, n):
            if (subset_mask >> j) & 1:
                result[i, j] *= -1
                result[j, i] *= -1
    return result


def core_orbit_pressure(matrix: np.ndarray, raw_temperature: float) -> float:
    """Evaluate log E_J Zbar_{A^J}(2t), the exact K_0 identity."""
    logs = []
    for subset_mask in range(1 << len(matrix)):
        values = energy_values(induced_flip(matrix, subset_mask))
        logs.append(log_mean_cosh(values, 2.0 * raw_temperature))
    peak = max(logs)
    return peak + math.log(math.fsum(math.exp(value - peak) for value in logs)) - (
        len(matrix) * math.log(2.0)
    )


def matrix_sha256(matrix: np.ndarray) -> str:
    return hashlib.sha256(matrix.tobytes()).hexdigest()


def audit(max_r: int, betas: list[float]) -> dict[str, object]:
    payload: dict[str, object] = {
        "schema": "quadratic-signing-universal-double-minimizers-v1",
        "status": (
            "exhaustive switching-gauge enumeration; floating-point evaluation "
            "of exact integer energy histograms"
        ),
        "normalization": {
            "energy": "H_A(x)=sum_{i<j} a_ij x_i x_j",
            "pressure": "phi_A(t)=log E_x cosh(t H_A(x))",
            "child_raw_temperature": "s=beta/sqrt(r)",
            "parent_raw_temperature": "t=beta/sqrt(2r)",
            "candidate_defect": "phi_K(t)-2 P_r(beta)",
        },
        "global_negation_policy": (
            "retained as a distinct class because K_+(A) can distinguish A and -A"
        ),
        "max_child_order": max_r,
        "betas": betas,
        "orders": [],
    }

    for n in range(1, max_r + 1):
        values_by_code, histogram_fibres = enumerate_root_energy_data(n)
        minimizer_codes_by_beta: dict[float, set[int]] = {}
        grid_records: list[dict[str, object]] = []

        for beta in betas:
            raw_temperature = beta / math.sqrt(n)
            scored = sorted(
                (
                    pressure_from_absolute_histogram(key, raw_temperature),
                    key,
                )
                for key in histogram_fibres
            )
            optimum = scored[0][0]
            minimizing_keys = [
                key for value, key in scored if abs(value - optimum) <= 1e-12
            ]
            selected_codes = {
                code
                for key in minimizing_keys
                for code in histogram_fibres[key]
            }
            minimizer_codes_by_beta[beta] = selected_codes
            next_gap = None
            for value, _key in scored:
                if value > optimum + 1e-12:
                    next_gap = value - optimum
                    break
            grid_records.append(
                {
                    "beta": beta,
                    "child_raw_temperature": raw_temperature,
                    "P_r_beta": optimum,
                    "gap_to_next_distinct_histogram": next_gap,
                    "minimizing_absolute_energy_histograms": [
                        histogram_record(key) for key in minimizing_keys
                    ],
                    "minimizing_root_gauge_code_count": len(selected_codes),
                }
            )

        selected_union = set().union(*minimizer_codes_by_beta.values())
        classes = permutation_classes(n, selected_union)
        free_bit_count = len(free_edges(n))
        complement_mask = (1 << free_bit_count) - 1
        class_records: list[dict[str, object]] = []

        for class_record in classes:
            code = int(class_record["canonical_root_code"])
            matrix = matrix_from_root_code(n, code)
            child_values = values_by_code[code]
            core = universal_double(matrix, matching=0)
            integral = universal_double(matrix, matching=1)
            core_values = energy_values(core)
            integral_values = energy_values(integral)
            beta_records = []
            maximum_identity_error = 0.0

            for beta in betas:
                active_codes = minimizer_codes_by_beta[beta]
                if not active_codes.intersection(
                    set(class_record["_selected_root_codes"])
                ):
                    continue
                child_t = beta / math.sqrt(n)
                parent_t = beta / math.sqrt(2 * n)
                child_pressure = log_mean_cosh(child_values, child_t)
                core_pressure = log_mean_cosh(core_values, parent_t)
                integral_pressure = log_mean_cosh(integral_values, parent_t)
                orbit_pressure = core_orbit_pressure(matrix, parent_t)
                identity_error = abs(core_pressure - orbit_pressure)
                if identity_error > 1e-10:
                    raise AssertionError((n, code, beta, identity_error))
                matching_bound = beta * math.sqrt(n / 2.0)
                if abs(integral_pressure - core_pressure) > matching_bound + 1e-10:
                    raise AssertionError(
                        (n, code, beta, integral_pressure - core_pressure)
                    )
                maximum_identity_error = max(maximum_identity_error, identity_error)
                beta_records.append(
                    {
                        "beta": beta,
                        "child_pressure": child_pressure,
                        "core_pressure": core_pressure,
                        "integral_plus_I_pressure": integral_pressure,
                        "core_orbit_identity_pressure": orbit_pressure,
                        "core_orbit_identity_abs_error": identity_error,
                        "core_defect_vs_2P_r": core_pressure - 2 * child_pressure,
                        "integral_plus_I_defect_vs_2P_r": (
                            integral_pressure - 2 * child_pressure
                        ),
                        "matching_pressure_difference": (
                            integral_pressure - core_pressure
                        ),
                        "matching_absolute_bound_beta_sqrt_r_over_2": (
                            matching_bound
                        ),
                    }
                )

            class_records.append(
                {
                    **{
                        key: value
                        for key, value in class_record.items()
                        if not key.startswith("_")
                    },
                    "global_negation_root_code": code ^ complement_mask,
                    "matrix": matrix.astype(int).tolist(),
                    "matrix_sha256": matrix_sha256(matrix),
                    "child_signed_energy_histogram_projective": signed_histogram(
                        child_values
                    ),
                    "child_cap": int(np.max(np.abs(child_values))),
                    "core_signed_energy_histogram_projective": signed_histogram(
                        core_values
                    ),
                    "core_cap": int(np.max(np.abs(core_values))),
                    "integral_plus_I_signed_energy_histogram_projective": (
                        signed_histogram(integral_values)
                    ),
                    "integral_plus_I_cap": int(np.max(np.abs(integral_values))),
                    "max_core_orbit_identity_abs_error": maximum_identity_error,
                    "beta_records": beta_records,
                }
            )

        # Add the best and worst construction values across nonunique exact
        # child classes.  This is the relevant selector range for a recurrence.
        for record in grid_records:
            beta = float(record["beta"])
            active = [
                beta_record
                for class_record in class_records
                for beta_record in class_record["beta_records"]
                if float(beta_record["beta"]) == beta
            ]
            record["switching_permutation_class_count"] = len(active)
            record["core_defect_range"] = [
                min(float(item["core_defect_vs_2P_r"]) for item in active),
                max(float(item["core_defect_vs_2P_r"]) for item in active),
            ]
            record["integral_plus_I_defect_range"] = [
                min(
                    float(item["integral_plus_I_defect_vs_2P_r"])
                    for item in active
                ),
                max(
                    float(item["integral_plus_I_defect_vs_2P_r"])
                    for item in active
                ),
            ]

        payload["orders"].append(
            {
                "r": n,
                "root_gauge_signings_enumerated": len(values_by_code),
                "distinct_absolute_energy_histograms": len(histogram_fibres),
                "all_absolute_energy_histogram_fibres": [
                    {
                        "histogram": histogram_record(key),
                        "root_gauge_code_count": len(histogram_fibres[key]),
                    }
                    for key in sorted(histogram_fibres)
                ],
                "grid": grid_records,
                "minimizer_classes": class_records,
            }
        )

    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-r", type=int, default=7)
    parser.add_argument("--betas", type=float, nargs="+", default=DEFAULT_BETAS)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "computations/results/universal_double_actual_minimizers_n7.json"
        ),
    )
    parser.add_argument("--quiet", action="store_true")
    args = parser.parse_args()
    result = audit(args.max_r, args.betas)
    rendered = json.dumps(result, indent=2, sort_keys=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered + "\n", encoding="utf-8")
    if not args.quiet:
        print(rendered)


if __name__ == "__main__":
    main()
