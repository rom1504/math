#!/usr/bin/env python3
"""Exact-cube audit of CC.14 and the IC.3 tilted average-influence path.

The child signings are the complete contracted-temperature minimizer classes
from ``actual_child_bridge_law_exact``.  Every bridge is enumerated.  This is
finite numerical evidence, not an interval certificate.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import sys
from pathlib import Path

import mpmath as mp
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import actual_child_bridge_law_exact as exact  # noqa: E402


def logsumexp(values: np.ndarray) -> float:
    maximum = float(np.max(values))
    return maximum + math.log(float(np.sum(np.exp(values - maximum))))


def logmeanexp(values: np.ndarray) -> float:
    return logsumexp(values) - math.log(len(values))


def row_marginal_logs(
    log_likelihood: np.ndarray, m: int, n: int, row: int
) -> np.ndarray:
    """Return log likelihood of one output row relative to its fair law."""

    lower_bits = row * n
    upper_bits = m * n - lower_bits - n
    view = log_likelihood.reshape(
        1 << upper_bits, 1 << n, 1 << lower_bits
    )
    result = np.empty(1 << n, dtype=np.float64)
    for code in range(1 << n):
        result[code] = logmeanexp(view[:, code, :].reshape(-1))
    return result


def projective_record(
    pressure: np.ndarray,
    m: int,
    n: int,
    lam: float,
    quadrature_order: int,
) -> dict:
    d = m * n
    log_p = pressure - logmeanexp(pressure)
    log_z = [row_marginal_logs(log_p, m, n, row) for row in range(m)]
    deltas = []
    context_oscillations = []
    for row in range(m):
        lower_bits = row * n
        upper_bits = d - lower_bits - n
        view = log_p.reshape(1 << upper_bits, 1 << n, 1 << lower_bits)
        score = view - log_z[row][None, :, None]
        context_oscillation = np.max(score, axis=1) - np.min(score, axis=1)
        context_oscillations.append(context_oscillation)
        deltas.append(float(np.max(context_oscillation)))

    masks = np.arange(1 << d, dtype=np.uint64)
    log_r = np.zeros(1 << d, dtype=np.float64)
    sum_log_z = np.zeros(1 << d, dtype=np.float64)
    row_log_probabilities = []
    row_mask = np.uint64((1 << n) - 1)
    for row in range(m):
        log_weights = -lam * log_z[row]
        log_prob = log_weights - logsumexp(log_weights)
        row_log_probabilities.append(log_prob)
        codes = ((masks >> np.uint64(row * n)) & row_mask).astype(np.int64)
        log_r += log_prob[codes]
        sum_log_z += log_z[row][codes]

    h = log_p - sum_log_z
    r_probability = np.exp(log_r)
    mean_h = float(np.dot(r_probability, h))
    canonical_j = logsumexp(log_r - lam * (h - mean_h))
    average_projective_squared = []
    conditional_variance_terms = []
    row_centered_differences = []
    for row in range(m):
        lower_bits = row * n
        upper_bits = d - lower_bits - n
        shape = (1 << upper_bits, 1 << n, 1 << lower_bits)
        context_probability = r_probability.reshape(shape).sum(axis=1)
        average_projective_squared.append(
            float(np.sum(context_probability * context_oscillations[row] ** 2))
        )
        h_view = h.reshape(shape)
        row_probability = np.exp(row_log_probabilities[row])
        conditional_mean = np.sum(
            h_view * row_probability[None, :, None], axis=1
        )
        conditional_second = np.sum(
            h_view * h_view * row_probability[None, :, None], axis=1
        )
        conditional_variance_terms.append(
            float(
                np.sum(
                    context_probability
                    * np.maximum(
                        conditional_second - conditional_mean * conditional_mean,
                        0.0,
                    )
                )
            )
        )
        row_centered_differences.append(
            (h_view - conditional_mean[:, None, :]).reshape(-1)
        )
    delta_squared = float(np.dot(deltas, deltas))
    theorem_bound = lam * lam * delta_squared / 8.0

    gauss_nodes, gauss_weights = np.polynomial.legendre.leggauss(
        quadrature_order
    )
    integration_nodes = 0.5 * lam * (gauss_nodes + 1.0)
    integration_weights = 0.5 * lam * gauss_weights
    profile_nodes = sorted(
        set([0.0, lam, *(float(value) for value in integration_nodes)])
    )
    path_profile = []
    path_by_s = {}
    for s in profile_nodes:
        log_q = log_r - s * h
        log_q -= logsumexp(log_q)
        q_probability = np.exp(log_q)
        q_mean_h = float(np.dot(q_probability, h))
        curvature = float(np.dot(q_probability, (h - q_mean_h) ** 2))
        divergence_qs_from_r = float(
            np.dot(q_probability, log_q - log_r)
        )
        row_marginal_drift = 0.0
        for row in range(m):
            lower_bits = row * n
            upper_bits = d - lower_bits - n
            shape = (1 << upper_bits, 1 << n, 1 << lower_bits)
            q_row = q_probability.reshape(shape).sum(axis=(0, 2))
            row_marginal_drift += float(
                np.dot(
                    q_row,
                    np.log(q_row) - row_log_probabilities[row],
                )
            )
        row_total_correlation = divergence_qs_from_r - row_marginal_drift
        if s == 0.0:
            influence = 0.5 * sum(conditional_variance_terms)
            entropic_influence = influence
        else:
            tau_sum = np.zeros_like(h)
            for difference in row_centered_differences:
                argument = s * difference
                tau_sum += np.expm1(argument) - argument
            influence = float(np.dot(q_probability, tau_sum)) / (s * s)
            entropic_kl_sum = 0.0
            for row in range(m):
                lower_bits = row * n
                upper_bits = d - lower_bits - n
                shape = (1 << upper_bits, 1 << n, 1 << lower_bits)
                h_view = h.reshape(shape)
                weighted = (
                    -s * h_view
                    + row_log_probabilities[row][None, :, None]
                )
                maximum = np.max(weighted, axis=1)
                log_normalizer = maximum + np.log(
                    np.sum(
                        np.exp(weighted - maximum[:, None, :]), axis=1
                    )
                )
                h_i_s = -log_normalizer / s
                q_view = q_probability.reshape(shape)
                entropic_kl_sum += float(
                    np.sum(q_view * (-s * (h_view - h_i_s[:, None, :])))
                )
            entropic_influence = entropic_kl_sum / (s * s)
        centered_cumulant = logsumexp(
            log_r - s * (h - mean_h)
        )
        item = {
            "s": s,
            "K_s": centered_cumulant,
            "curvature_Var_qs_h": curvature,
            "curvature_per_parent_vertex": curvature / (m + n),
            "tilted_average_influence_A_s": influence,
            "tilted_average_influence_per_parent_vertex": influence
            / (m + n),
            "conditional_entropic_influence_E_s": entropic_influence,
            "conditional_entropic_influence_per_parent_vertex": (
                entropic_influence / (m + n)
            ),
            "D_qs_from_canonical_row_product": divergence_qs_from_r,
            "row_total_correlation_qs": row_total_correlation,
            "row_marginal_drift_from_canonical_factors": row_marginal_drift,
            "row_total_correlation_fraction_of_D_qs_r": (
                row_total_correlation / divergence_qs_from_r
                if divergence_qs_from_r > 1e-14
                else None
            ),
            "row_marginal_drift_fraction_of_D_qs_r": (
                row_marginal_drift / divergence_qs_from_r
                if divergence_qs_from_r > 1e-14
                else None
            ),
        }
        path_profile.append(item)
        path_by_s[s] = item

    curvature_quadrature = 0.0
    influence_quadrature = 0.0
    entropic_influence_quadrature = 0.0
    total_correlation_quadrature = 0.0
    row_marginal_drift_quadrature = 0.0
    for s, weight in zip(integration_nodes, integration_weights):
        item = path_by_s[float(s)]
        curvature_quadrature += (
            weight * (lam - s) * item["curvature_Var_qs_h"]
        )
        influence_quadrature += (
            weight * item["tilted_average_influence_A_s"]
        )
        entropic_influence_quadrature += (
            weight * item["conditional_entropic_influence_E_s"]
        )
        total_correlation_quadrature += (
            weight * item["row_total_correlation_qs"] / (s * s)
        )
        row_marginal_drift_quadrature += (
            weight
            * item["row_marginal_drift_from_canonical_factors"]
            / (s * s)
        )
    ic15_rhs = lam * influence_quadrature
    ic23_rhs = lam * entropic_influence_quadrature
    es28_total_correlation = lam * total_correlation_quadrature
    es28_row_marginal_drift = lam * row_marginal_drift_quadrature
    es28_total = es28_total_correlation + es28_row_marginal_drift
    marginal_spread = max(
        float(np.max(np.abs(log_z[row] - log_z[0])))
        for row in range(m)
    )
    return {
        "row_projective_diameters": deltas,
        "projective_Delta_squared": delta_squared,
        "projective_Delta_squared_per_parent_vertex": delta_squared / (m + n),
        "row_average_projective_squared": average_projective_squared,
        "sum_average_projective_squared": sum(average_projective_squared),
        "sum_average_projective_squared_per_parent_vertex": sum(
            average_projective_squared
        )
        / (m + n),
        "row_conditional_variance_terms": conditional_variance_terms,
        "Efron_Stein_variance_upper_bound": sum(conditional_variance_terms),
        "actual_variance_h_under_canonical_product": float(
            np.dot(r_probability, (h - mean_h) ** 2)
        ),
        "interaction_path_quadrature_order": quadrature_order,
        "interaction_path_profile": path_profile,
        "maximum_tilted_average_influence": max(
            item["tilted_average_influence_A_s"] for item in path_profile
        ),
        "maximum_tilted_average_influence_per_parent_vertex": max(
            item["tilted_average_influence_per_parent_vertex"]
            for item in path_profile
        ),
        "maximum_conditional_entropic_influence": max(
            item["conditional_entropic_influence_E_s"]
            for item in path_profile
        ),
        "maximum_conditional_entropic_influence_per_parent_vertex": max(
            item["conditional_entropic_influence_per_parent_vertex"]
            for item in path_profile
        ),
        "IC7_curvature_quadrature": curvature_quadrature,
        "IC7_curvature_quadrature_to_exact_J_ratio": (
            curvature_quadrature / canonical_j
            if canonical_j > 1e-14
            else None
        ),
        "IC15_rhs_quadrature": ic15_rhs,
        "IC15_rhs_quadrature_to_exact_J_ratio": (
            ic15_rhs / canonical_j if canonical_j > 1e-14 else None
        ),
        "IC23_rhs_quadrature": ic23_rhs,
        "IC23_rhs_quadrature_to_exact_J_ratio": (
            ic23_rhs / canonical_j if canonical_j > 1e-14 else None
        ),
        "ES28_total_correlation_contribution_quadrature": (
            es28_total_correlation
        ),
        "ES28_row_marginal_drift_contribution_quadrature": (
            es28_row_marginal_drift
        ),
        "ES28_decomposition_quadrature": es28_total,
        "ES28_decomposition_quadrature_to_exact_J_ratio": (
            es28_total / canonical_j if canonical_j > 1e-14 else None
        ),
        "ES28_total_correlation_fraction_of_decomposition": (
            es28_total_correlation / es28_total
            if es28_total > 1e-14
            else None
        ),
        "ES28_row_marginal_drift_fraction_of_decomposition": (
            es28_row_marginal_drift / es28_total
            if es28_total > 1e-14
            else None
        ),
        "canonical_J": canonical_j,
        "canonical_J_per_parent_vertex": canonical_j / (m + n),
        "CC2_upper_bound": theorem_bound,
        "CC2_upper_bound_per_parent_vertex": theorem_bound / (m + n),
        "bound_to_J_ratio": theorem_bound / canonical_j
        if canonical_j > 1e-14
        else None,
        "maximum_row_marginal_log_spread": marginal_spread,
        "row_product_log_probability_normalization_error": abs(
            logsumexp(log_r)
        ),
    }


def children(order: int, beta: float, total_n: int) -> list[dict]:
    space = exact.build_signing_space(order)
    classes, selector = exact.thermal_minimizer_classes(
        space, beta, total_n
    )
    return [
        {
            "matrix": np.asarray(item["representative_matrix"], dtype=np.int16),
            "sha256": item["representative_sha256"],
            "selector": selector,
        }
        for item in classes
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orders", nargs="+", type=int, default=list(range(4, 10)))
    parser.add_argument("--betas", nargs="+", type=float, default=[1.0, 2.0, 4.0])
    parser.add_argument("--lambda-value", type=float, default=1.0)
    parser.add_argument("--quadrature-order", type=int, default=12)
    parser.add_argument(
        "--mp-dps",
        type=int,
        default=80,
        help=(
            "decimal precision for the optimizing-child histogram comparison; "
            "must exceed the 20 guard digits used by the selector"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=(
            ROOT
            / "computations/results/actual_child_projective_synchronization.json"
        ),
    )
    args = parser.parse_args()
    if args.mp_dps <= 20:
        parser.error("--mp-dps must be greater than 20")
    mp.mp.dps = args.mp_dps

    records = []
    for total_n in args.orders:
        m = total_n // 2
        n = total_n - m
        for beta in args.betas:
            left_children = children(m, beta, total_n)
            right_children = children(n, beta, total_n)
            for left_index, right_index, epsilon in itertools.product(
                range(len(left_children)),
                range(len(right_children)),
                (-1, 1),
            ):
                left = left_children[left_index]
                right = right_children[right_index]
                pressure, pressure_audit = exact.bridge_pressures(
                    left["matrix"],
                    right["matrix"],
                    beta,
                    total_n,
                    epsilon,
                )
                records.append(
                    {
                        "N": total_n,
                        "split": [m, n],
                        "beta": beta,
                        "lambda": args.lambda_value,
                        "epsilon": epsilon,
                        "left_class": left_index,
                        "right_class": right_index,
                        "left_sha256": left["sha256"],
                        "right_sha256": right["sha256"],
                        "pressure_audit": pressure_audit,
                        **projective_record(
                            pressure,
                            m,
                            n,
                            args.lambda_value,
                            args.quadrature_order,
                        ),
                    }
                )
                print(
                    f"N={total_n} beta={beta:g} "
                    f"classes={left_index},{right_index} eps={epsilon:+d} "
                    f"Delta2/N={records[-1]['projective_Delta_squared_per_parent_vertex']:.6g} "
                    f"J/N={records[-1]['canonical_J_per_parent_vertex']:.6g}",
                    flush=True,
                )

    payload = {
        "schema": "actual-child-projective-synchronization-and-influence-v3",
        "classification": "complete finite bridge enumeration; numerical, not interval-certified",
        "parameters": {
            "orders": args.orders,
            "betas": args.betas,
            "lambda": args.lambda_value,
            "mp_dps": args.mp_dps,
            "interaction_path_gauss_legendre_order": args.quadrature_order,
            "balanced_split_only": True,
            "all_contracted_temperature_minimizer_classes": True,
            "all_relative_orientations": True,
        },
        "records": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
