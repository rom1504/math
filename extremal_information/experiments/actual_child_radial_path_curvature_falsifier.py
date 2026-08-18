#!/usr/bin/env python3
"""Actual-minimizer falsifier for radial FC data versus the ES path.

The exact part computes the signed ground-state tangent matrices of the two
certified order-eight minimizer classes.  The numerical part attaches the
unique order-two child at raw temperature/channel amplitude three and
evaluates the complete 2 x 8 bridge cube.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from fractions import Fraction
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))

import actual_child_bridge_law_exact as exact  # noqa: E402
import actual_child_projective_synchronization as path  # noqa: E402


def energy_histogram(energies: np.ndarray) -> dict[str, int]:
    values, counts = np.unique(energies, return_counts=True)
    return {str(int(value)): int(count) for value, count in zip(values, counts)}


def signed_tangent(matrix: np.ndarray, raw_t: float) -> np.ndarray:
    spins = exact.projective_spins(len(matrix)).astype(np.float64)
    energies = exact.energies_for_matrix(matrix, spins.astype(np.int16))
    denominator = float(np.mean(np.cosh(raw_t * energies)))
    return (
        np.einsum(
            "b,bi,bj->ij",
            np.sinh(raw_t * energies),
            spins,
            spins,
        )
        / len(spins)
        / denominator
    )


def exact_ground_tangent(matrix: np.ndarray) -> tuple[np.ndarray, int, int]:
    spins = exact.projective_spins(len(matrix)).astype(np.int16)
    energies = exact.energies_for_matrix(matrix, spins)
    cap = int(np.max(np.abs(energies)))
    active = np.abs(energies) == cap
    signed_sum = np.einsum(
        "b,bi,bj->ij",
        np.sign(energies[active]).astype(np.int64),
        spins[active].astype(np.int64),
        spins[active].astype(np.int64),
        dtype=np.int64,
    )
    return signed_sum, int(np.count_nonzero(active)), cap


def physical_tropical_certificate(
    left: np.ndarray, right: np.ndarray, epsilon: int
) -> dict:
    """Exact lambda=1 slope of J/t when internal and bridge t tend to infinity."""

    m, n = len(left), len(right)
    if m != 2:
        raise ValueError("the exact rational certificate is specialized to two rows")
    x = exact.projective_spins(m).astype(np.int16)
    y = exact.projective_spins(n).astype(np.int16)
    ex = exact.energies_for_matrix(left, x)
    ey = exact.energies_for_matrix(right, y)
    internal = np.abs((ex[:, None] + epsilon * ey[None, :]).reshape(-1))
    rank_one = np.asarray(
        [(xx[:, None] * yy[None, :]).reshape(-1) for xx in x for yy in y],
        dtype=np.int16,
    )
    d = m * n
    masks = np.arange(1 << d, dtype=np.uint64)
    bits = (
        (masks[:, None] >> np.arange(d, dtype=np.uint64)) & 1
    ).astype(np.int16)
    bridges = 1 - 2 * bits
    cross = np.abs(bridges @ rank_one.T)
    rates = cross + internal[None, :]
    bridge_rate = np.max(rates, axis=1)
    active_count = np.sum(rates == bridge_rate[:, None], axis=1)
    # Every maximizing term in this witness has both cosh arguments nonzero,
    # so the common leading coefficient is active_count/4.
    active = rates == bridge_rate[:, None]
    if np.any(active & (internal[None, :] == 0)) or np.any(
        active & (cross == 0)
    ):
        raise AssertionError("unexpected unequal leading cosh coefficient")

    row_size = 1 << n
    rate_table = bridge_rate.reshape(row_size, row_size)
    count_table = active_count.reshape(row_size, row_size)
    lower_rate = np.max(rate_table, axis=0)
    upper_rate = np.max(rate_table, axis=1)
    lower_support = np.flatnonzero(lower_rate == np.min(lower_rate))
    upper_support = np.flatnonzero(upper_rate == np.min(upper_rate))
    lower_coefficient = np.asarray(
        [
            np.sum(count_table[rate_table[:, j] == lower_rate[j], j])
            for j in range(row_size)
        ],
        dtype=np.int64,
    )
    upper_coefficient = np.asarray(
        [
            np.sum(count_table[i, rate_table[i, :] == upper_rate[i]])
            for i in range(row_size)
        ],
        dtype=np.int64,
    )
    lower_weights = [Fraction(1, int(lower_coefficient[j])) for j in lower_support]
    upper_weights = [Fraction(1, int(upper_coefficient[i])) for i in upper_support]
    lower_normalizer = sum(lower_weights, Fraction(0))
    upper_normalizer = sum(upper_weights, Fraction(0))
    expected_rate = Fraction(0)
    for i, wi in zip(upper_support, upper_weights):
        for j, wj in zip(lower_support, lower_weights):
            expected_rate += (
                wi
                * wj
                * int(rate_table[i, j])
                / (upper_normalizer * lower_normalizer)
            )
    minimum_rate = int(np.min(rate_table))
    slope = expected_rate - minimum_rate
    return {
        "orientation": epsilon,
        "global_minimum_rate": minimum_rate,
        "global_maximum_rate": int(np.max(rate_table)),
        "lower_row_minimax_rate": int(np.min(lower_rate)),
        "upper_row_minimax_rate": int(np.min(upper_rate)),
        "lower_row_tropical_support_size": len(lower_support),
        "upper_row_tropical_support_size": len(upper_support),
        "lower_row_support_active_coefficients": [
            [int(j), int(lower_coefficient[j])] for j in lower_support
        ],
        "upper_row_support_active_coefficients": [
            [int(i), int(upper_coefficient[i])] for i in upper_support
        ],
        "canonical_r_expected_parent_rate": {
            "numerator": expected_rate.numerator,
            "denominator": expected_rate.denominator,
        },
        "exact_limit_J_over_t_at_lambda_1": {
            "numerator": slope.numerator,
            "denominator": slope.denominator,
            "float": float(slope),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-t", type=float, default=3.0)
    parser.add_argument("--lambda-value", type=float, default=1.0)
    parser.add_argument("--quadrature-order", type=int, default=12)
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            ROOT
            / "computations/results/actual_child_radial_path_curvature_falsifier.json"
        ),
    )
    args = parser.parse_args()

    classification_path = ROOT / "computations/results/m8_minimizer_orbits.json"
    classification = json.loads(classification_path.read_text(encoding="utf-8"))
    if classification["order"] != 8 or classification["target_cap"] != 10:
        raise AssertionError("unexpected order-eight minimizer certificate")
    if len(classification["classes"]) != 2:
        raise AssertionError("the falsifier requires the two certified classes")

    right_children = [
        np.asarray(item["representative_matrix"], dtype=np.int16)
        for item in classification["classes"]
    ]
    left_child = np.asarray([[0, 1], [1, 0]], dtype=np.int16)
    total_n = 10
    beta = args.raw_t * math.sqrt(total_n)

    records = []
    histograms = []
    for class_index, right in enumerate(right_children):
        spins = exact.projective_spins(8).astype(np.int16)
        energies = exact.energies_for_matrix(right, spins)
        histogram = energy_histogram(energies)
        histograms.append(histogram)

        ground_sum, ground_count, cap = exact_ground_tangent(right)
        ground_trace_numerator = int(np.sum(ground_sum * ground_sum))
        ground_trace_denominator = ground_count * ground_count
        if ground_trace_numerator % ground_trace_denominator:
            raise AssertionError("expected integral limiting trace")
        ground_trace = ground_trace_numerator // ground_trace_denominator

        tangent = signed_tangent(right, args.raw_t)
        tangent_trace = float(np.sum(tangent * tangent))
        gamma_squared = math.tanh(args.raw_t) ** 2 * tangent_trace
        leading_j_coefficient = 0.5 * args.lambda_value**2 * gamma_squared

        physical_records = []
        for epsilon in (-1, 1):
            pressure, pressure_audit = exact.bridge_pressures(
                left_child, right, beta, total_n, epsilon
            )
            response = path.projective_record(
                pressure,
                2,
                8,
                args.lambda_value,
                args.quadrature_order,
            )
            physical_records.append(
                {
                    "epsilon": epsilon,
                    "pressure_audit": pressure_audit,
                    "canonical_J": response["canonical_J"],
                    "canonical_J_per_parent_vertex": response[
                        "canonical_J_per_parent_vertex"
                    ],
                    "ES28_total_correlation_fraction": response[
                        "ES28_total_correlation_fraction_of_decomposition"
                    ],
                    "ES28_row_marginal_drift_fraction": response[
                        "ES28_row_marginal_drift_fraction_of_decomposition"
                    ],
                    "IC23_to_J_ratio": response[
                        "IC23_rhs_quadrature_to_exact_J_ratio"
                    ],
                    "ES28_to_J_ratio": response[
                        "ES28_decomposition_quadrature_to_exact_J_ratio"
                    ],
                }
            )

        records.append(
            {
                "class_index": class_index,
                "matrix_sha256": classification["classes"][class_index][
                    "representative_matrix_sha256"
                ],
                "energy_histogram": histogram,
                "cap": cap,
                "ground_state_count": ground_count,
                "exact_signed_ground_tangent_numerator": ground_sum.tolist(),
                "exact_signed_ground_tangent_denominator": ground_count,
                "exact_limiting_tangent_frobenius_squared": ground_trace,
                "finite_t_tangent_frobenius_squared": tangent_trace,
                "finite_t_cross_row_Gamma_frobenius_squared": gamma_squared,
                "low_channel_rho4_leading_J_coefficient": leading_j_coefficient,
                "physical_tropical_certificate": physical_tropical_certificate(
                    left_child, right, -1
                ),
                "physical_channel_records": physical_records,
            }
        )

    if histograms[0] != histograms[1]:
        raise AssertionError("the two radial pressure histograms must agree")
    if [item["exact_limiting_tangent_frobenius_squared"] for item in records] != [
        14,
        10,
    ]:
        raise AssertionError("unexpected exact ground tangent traces")

    payload = {
        "schema": "actual-child-radial-path-curvature-falsifier-v1",
        "classification": (
            "exact finite ground-state arithmetic and certified minimizer input; "
            "complete floating bridge enumeration at raw t=3 is numerical, not interval-certified"
        ),
        "parameters": {
            "internal_raw_t": args.raw_t,
            "physical_channel_u": args.raw_t,
            "beta_for_parent_normalization": beta,
            "lambda": args.lambda_value,
            "total_parent_order": total_n,
            "split": [2, 8],
            "quadrature_order": args.quadrature_order,
        },
        "minimizer_certificate": {
            "source": str(classification_path.relative_to(ROOT)),
            "order": classification["order"],
            "target_cap": classification["target_cap"],
            "class_count": classification[
                "signed_permutation_and_global_sign_class_count"
            ],
            "common_histogram": histograms[0],
            "large_temperature_fact": (
                "FC.22 proves both classes minimize pressure for every raw t>=3"
            ),
        },
        "records": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote {args.output}")
    for record in records:
        physical = record["physical_channel_records"][0]
        print(
            f"class={record['class_index']} "
            f"ground_trace={record['exact_limiting_tangent_frobenius_squared']} "
            f"J/N={physical['canonical_J_per_parent_vertex']:.9f} "
            f"TCshare={physical['ES28_total_correlation_fraction']:.9f}",
            flush=True,
        )


if __name__ == "__main__":
    main()
