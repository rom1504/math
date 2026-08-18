#!/usr/bin/env python3
"""Finite actual-law check of the signed sector-mode product certificate.

Every child is selected at 80-digit precision before any bridge calculation.
The bridge cubes and the declared one-parameter product trial are then
evaluated completely.  Transcendental and eigensystem calculations are
floating point, so the output is numerical evidence, not an interval
certificate or an asymptotic theorem.
"""

from __future__ import annotations

import argparse
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
import actual_child_projective_synchronization as path  # noqa: E402
import actual_child_row_product_shadow as shadow  # noqa: E402


def sector_bias(matrix: np.ndarray, raw_t: float) -> float:
    spins = exact.projective_spins(len(matrix))
    energy = exact.energies_for_matrix(matrix, spins)
    positive = float(np.mean(np.exp(raw_t * energy)))
    negative = float(np.mean(np.exp(-raw_t * energy)))
    return 0.5 * math.log(positive / negative)


def unique_minimizer(order: int, beta: float, total_n: int) -> tuple[np.ndarray, dict]:
    classes, selector = exact.thermal_minimizer_classes(
        exact.build_signing_space(order), format(beta, ".12g"), total_n
    )
    if len(classes) != 1:
        raise AssertionError(
            f"expected one minimizing class at order={order}, N={total_n}, "
            f"beta={beta:g}; found {len(classes)}"
        )
    return (
        np.asarray(classes[0]["representative_matrix"], dtype=np.int8),
        selector,
    )


def balanced_children(
    total_n: int, beta: float
) -> tuple[np.ndarray, np.ndarray, int, float, float, list[dict]]:
    first_order = total_n // 2
    second_order = total_n - first_order
    first, first_selector = unique_minimizer(first_order, beta, total_n)
    second, second_selector = unique_minimizer(second_order, beta, total_n)
    raw_t = beta / math.sqrt(total_n)
    first_bias = sector_bias(first, raw_t)
    second_bias = sector_bias(second, raw_t)
    if abs(first_bias) > abs(second_bias) + 1e-13:
        first, second = second, first
        first_bias, second_bias = second_bias, first_bias
        first_selector, second_selector = second_selector, first_selector
    epsilon = -1 if first_bias * second_bias >= 0 else 1
    return (
        first,
        second,
        epsilon,
        first_bias,
        second_bias,
        [first_selector, second_selector],
    )


def signed_quadratic_matrix(
    left: np.ndarray,
    right: np.ndarray,
    raw_t: float,
    epsilon: int,
) -> np.ndarray:
    """Return the block matrix M in H_2(B)=B^T M B/2."""

    m, n = len(left), len(right)
    x = exact.projective_spins(m).astype(np.float64)
    y = exact.projective_spins(n).astype(np.float64)
    ex = exact.energies_for_matrix(left, x.astype(np.int16))
    ey = exact.energies_for_matrix(right, y.astype(np.int16))
    weight = np.cosh(raw_t * (ex[:, None] + epsilon * ey[None, :]))
    weight /= np.sum(weight)
    matrix = np.zeros((m * n, m * n), dtype=np.float64)
    for i in range(m):
        for k in range(i + 1, m):
            correlation = np.einsum(
                "a,b,ab,bj,bl->jl",
                x[:, i] * x[:, k],
                np.ones(len(y)),
                weight,
                y,
                y,
                optimize=True,
            )
            matrix[i * n : (i + 1) * n, k * n : (k + 1) * n] = correlation
            matrix[k * n : (k + 1) * n, i * n : (i + 1) * n] = correlation.T
    return matrix


def canonical_principal_vector(
    eigenvalues: np.ndarray, eigenvectors: np.ndarray
) -> tuple[np.ndarray, int, float | None]:
    """Choose a basis-invariant vector from the bottom eigenspace.

    Project the first coordinate vector having nonzero projection into the
    bottom eigenspace, then normalize.  This removes the arbitrary basis
    choice returned for a multiple eigenvalue.
    """

    scale = max(1.0, abs(float(eigenvalues[0])))
    eigenspace = np.flatnonzero(
        np.abs(eigenvalues - eigenvalues[0]) <= 1e-9 * scale
    )
    basis = eigenvectors[:, eigenspace]
    projector = basis @ basis.T
    vector = None
    for coordinate in range(len(eigenvalues)):
        candidate = projector[:, coordinate]
        norm = float(np.linalg.norm(candidate))
        if norm > 1e-9:
            vector = candidate / norm
            break
    if vector is None:
        raise AssertionError("empty principal eigenspace")
    first_nonzero = np.flatnonzero(np.abs(vector) > 1e-12)
    if len(first_nonzero) and vector[first_nonzero[0]] < 0:
        vector = -vector
    next_index = int(eigenspace[-1]) + 1
    gap = (
        float(eigenvalues[next_index] - eigenvalues[0])
        if next_index < len(eigenvalues)
        else None
    )
    return vector, len(eigenspace), gap


def row_words(n: int) -> np.ndarray:
    masks = np.arange(1 << n, dtype=np.uint64)
    bits = ((masks[:, None] >> np.arange(n, dtype=np.uint64)) & 1).astype(
        np.int8
    )
    return (1 - 2 * bits).astype(np.float64)


def record(total_n: int, beta: float, lam: float) -> dict:
    raw_t = beta / math.sqrt(total_n)
    left, right, epsilon, left_bias, right_bias, selectors = balanced_children(
        total_n, beta
    )
    m, n = len(left), len(right)
    matrix = signed_quadratic_matrix(left, right, raw_t, epsilon)
    eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    mode, multiplicity, eigengap = canonical_principal_vector(
        eigenvalues, eigenvectors
    )

    pressure, pressure_audit = exact.bridge_pressures(
        left, right, beta, total_n, epsilon
    )
    log_p = pressure - path.logmeanexp(pressure)
    row_logs = [
        path.row_marginal_logs(log_p, m, n, row) for row in range(m)
    ]
    row_spread = max(float(np.max(np.abs(item - row_logs[0]))) for item in row_logs)
    if row_spread > 2e-12:
        raise AssertionError(f"canonical rows are not iid: {row_spread}")
    log_z = row_logs[0]
    row_probability = shadow.softmax(-lam * log_z)

    words = row_words(n)
    features: list[np.ndarray | None] = []
    response_vectors = []
    active_rows = []
    feature_fairness = []
    tie_counts = []
    for row in range(m):
        block = mode[row * n : (row + 1) * n]
        if np.linalg.norm(block) <= 1e-9:
            features.append(None)
            response_vectors.append(np.zeros(n, dtype=np.float64))
            continue
        active_rows.append(row)
        score = words @ block
        feature = np.sign(score)
        ties = np.flatnonzero(np.abs(score) <= 1e-12)
        feature[ties] = words[ties, 0]
        fairness = float(np.dot(row_probability, feature))
        if abs(fairness) > 2e-12:
            raise AssertionError(f"odd feature is not fair: {fairness}")
        features.append(feature)
        response_vectors.append((row_probability * feature) @ words)
        feature_fairness.append(fairness)
        tie_counts.append(len(ties))

    active_count = len(active_rows)
    if active_count == 0:
        raise AssertionError("principal mode has no active row")
    response = np.concatenate(response_vectors)
    signed_energy = -lam * raw_t**2 * float(response @ matrix @ response)
    rounded_ratio = signed_energy / active_count

    masks = np.arange(len(log_p), dtype=np.uint64)
    row_mask = np.uint64((1 << n) - 1)
    sum_log_z = np.zeros(len(log_p), dtype=np.float64)
    log_r = np.zeros(len(log_p), dtype=np.float64)
    codes = []
    for row in range(m):
        code = ((masks >> np.uint64(row * n)) & row_mask).astype(np.int64)
        codes.append(code)
        sum_log_z += log_z[code]
        log_r += np.log(row_probability[code])
    interaction = log_p - sum_log_z
    r_probability = np.exp(log_r)
    mean_r_interaction = float(np.dot(r_probability, interaction))
    centered = -lam * (interaction - mean_r_interaction)
    canonical_j = path.logsumexp(log_r + centered)

    trial = None
    if rounded_ratio > 1.0 + 1e-10:
        delta = min(rounded_ratio - 1.0, 1.0)
        amplitude = math.sqrt(3.0 * delta / (4.0 * (1.0 + delta)))
        log_trial = log_r.copy()
        for row in active_rows:
            feature = features[row]
            if feature is None:
                raise AssertionError("active row lost its feature")
            log_trial += amplitude * feature[codes[row]] - math.log(
                math.cosh(amplitude)
            )
        trial_probability = np.exp(log_trial)
        normalization_error = abs(float(np.sum(trial_probability)) - 1.0)
        if normalization_error > 2e-11:
            raise AssertionError(f"trial does not normalize: {normalization_error}")
        mean_trial_interaction = float(np.dot(trial_probability, interaction))
        binary_entropy_cost = active_count * (
            amplitude * math.tanh(amplitude) - math.log(math.cosh(amplitude))
        )
        actual_gain = (
            lam * (mean_r_interaction - mean_trial_interaction)
            - binary_entropy_cost
        )
        quadratic_gain = (
            0.5 * signed_energy * math.tanh(amplitude) ** 2
            - binary_entropy_cost
        )
        theorem_lower = (
            3.0 * delta**2 * active_count / (16.0 * (1.0 + delta))
        )
        trial = {
            "delta": delta,
            "amplitude": amplitude,
            "binary_entropy_cost": binary_entropy_cost,
            "quadratic_directional_gain": quadratic_gain,
            "SM25_quadratic_lower_bound": theorem_lower,
            "actual_directional_gain": actual_gain,
            "physical_remainder_contribution_to_gain": actual_gain
            - quadratic_gain,
            "trial_probability_normalization_error": normalization_error,
        }

    selector_gaps = [item["mp_gap_to_next_histogram"] for item in selectors]
    return {
        "N": total_n,
        "split": [m, n],
        "beta": beta,
        "lambda": lam,
        "raw_t": raw_t,
        "orientation": epsilon,
        "left_sha256": exact.matrix_sha(left),
        "right_sha256": exact.matrix_sha(right),
        "left_sector_bias": left_bias,
        "right_sector_bias": right_bias,
        "selector_mp_dps": mp.mp.dps,
        "selector_gaps_to_next_histogram": selector_gaps,
        "bridge_pressure_audit": pressure_audit,
        "canonical_row_log_spread": row_spread,
        "quadratic_K": 0.5 * float(np.sum(matrix * matrix)),
        "quadratic_lambda_min": float(eigenvalues[0]),
        "quadratic_lambda_max": float(eigenvalues[-1]),
        "principal_eigenspace_multiplicity": multiplicity,
        "principal_eigenspace_gap": eigengap,
        "active_mode_rows": active_rows,
        "active_mode_row_count": active_count,
        "feature_tie_counts": tie_counts,
        "maximum_feature_fairness_error": max(
            [abs(value) for value in feature_fairness], default=0.0
        ),
        "rounded_negative_energy_per_active_row": rounded_ratio,
        "canonical_J": canonical_j,
        "trial": trial,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orders", type=int, nargs="+", default=list(range(4, 10)))
    parser.add_argument("--betas", type=float, nargs="+", default=[1.0, 2.0, 4.0])
    parser.add_argument("--lambda-value", type=float, default=1.0)
    parser.add_argument("--mp-dps", type=int, default=80)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT
        / "computations/results/actual_child_signed_sector_mode_falsifier.json",
    )
    args = parser.parse_args()
    if args.mp_dps <= 20:
        parser.error("--mp-dps must be greater than 20")
    mp.mp.dps = args.mp_dps

    records = []
    for total_n in args.orders:
        for beta in args.betas:
            item = record(total_n, beta, args.lambda_value)
            records.append(item)
            trial = item["trial"]
            trial_text = (
                "no unstable rounded trial"
                if trial is None
                else (
                    f"quad={trial['quadratic_directional_gain']:.6g} "
                    f"actual={trial['actual_directional_gain']:.6g}"
                )
            )
            print(
                f"N={total_n} beta={beta:g} eps={item['orientation']:+d} "
                f"ratio={item['rounded_negative_energy_per_active_row']:.6g} "
                f"{trial_text}",
                flush=True,
            )

    payload = {
        "schema": "actual-child-signed-sector-mode-falsifier-v1",
        "classification": (
            "80-digit optimizing-child selection; complete finite bridge "
            "enumeration; floating transcendental/eigensystem evaluation; "
            "numerical evidence only"
        ),
        "parameters": {
            "orders": args.orders,
            "betas": args.betas,
            "lambda": args.lambda_value,
            "mp_dps": args.mp_dps,
            "balanced_orientation_and_row_direction": True,
            "canonical_principal_eigenspace_projection": True,
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
