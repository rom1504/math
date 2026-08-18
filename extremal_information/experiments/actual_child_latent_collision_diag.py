#!/usr/bin/env python3
"""Finite complete-cube diagnostic for actual-child latent posterior collision.

The contracted-temperature child prior is kept fixed while the bridge-channel
temperatures t and 2t are compared, as required by the exact posterior-
collision identity.  All displayed collision statistics are floating finite
evaluations; they are not asymptotic claims.
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
EXP = ROOT / "extremal_information" / "experiments"
sys.path.insert(0, str(EXP))

import actual_child_bridge_law_exact as exact  # noqa: E402
import actual_child_escort_low_degree_falsifier as low_degree  # noqa: E402


BETA = 4.0
LAMBDA = 1.0
PLANS = ((6, 3, 3), (8, 4, 4), (10, 3, 7))
ORIENTATIONS = (-1, 1)
QUANTILES = (0.01, 0.10, 0.50, 0.90, 0.99, 0.999)
CONTRIBUTION_LEVELS = (0.50, 0.90, 0.99)
NORMALIZED_LOG_THRESHOLDS = (0.10, 0.25, 0.50, 0.75, 1.00)


def logsumexp_weighted(log_values: np.ndarray, weights: np.ndarray) -> float:
    positive = weights > 0
    lv = log_values[positive]
    ww = weights[positive]
    maximum = float(np.max(lv))
    return maximum + math.log(float(np.dot(ww, np.exp(lv - maximum))))


def weighted_quantiles(
    values: np.ndarray, weights: np.ndarray, probabilities=QUANTILES
) -> dict[str, float]:
    order = np.argsort(values, kind="stable")
    cumulative = np.cumsum(weights[order], dtype=np.float64)
    cumulative /= cumulative[-1]
    answer = {}
    for probability in probabilities:
        index = min(
            len(order) - 1,
            int(np.searchsorted(cumulative, probability, side="left")),
        )
        answer[format(probability, ".12g")] = float(values[order[index]])
    return answer


def collision_summary(
    log_values: np.ndarray, weights: np.ndarray, total_order: int
) -> dict:
    weights = np.asarray(weights, dtype=np.float64)
    weights /= float(np.sum(weights))
    log_values = np.asarray(log_values, dtype=np.float64)
    log_mean = logsumexp_weighted(log_values, weights)

    descending = np.argsort(-log_values, kind="stable")
    log_tilt = (
        np.log(weights[descending], where=weights[descending] > 0,
               out=np.full(len(weights), -np.inf))
        + log_values[descending]
        - log_mean
    )
    contribution = np.exp(log_tilt)
    contribution_cumulative = np.cumsum(contribution)
    base_cumulative = np.cumsum(weights[descending])
    masses = {}
    for level in CONTRIBUTION_LEVELS:
        index = min(
            len(descending) - 1,
            int(np.searchsorted(contribution_cumulative, level, side="left")),
        )
        masses[format(level, ".12g")] = float(base_cumulative[index])

    # D(size-biased collision law || base law) is a coordinate-free measure
    # of how rare the words carrying the mean collision factor are.
    tilted = contribution / float(np.sum(contribution))
    divergence = float(np.dot(tilted, log_values[descending] - log_mean))

    return {
        "log_mean_collision": log_mean,
        "log_mean_collision_over_N": log_mean / total_order,
        "base_weighted_mean_log_collision_over_N": float(
            np.dot(weights, log_values) / total_order
        ),
        "base_weighted_log_collision_quantiles_over_N": {
            key: value / total_order
            for key, value in weighted_quantiles(log_values, weights).items()
        },
        "minimum_log_collision_over_N": float(np.min(log_values)) / total_order,
        "maximum_log_collision_over_N": float(np.max(log_values)) / total_order,
        "base_mass_log_collision_over_N_above": {
            format(threshold, ".12g"): float(
                np.sum(weights[log_values > threshold * total_order])
            )
            for threshold in NORMALIZED_LOG_THRESHOLDS
        },
        "base_mass_needed_for_collision_mean_contribution": masses,
        "collision_size_bias_divergence_from_base_nats": divergence,
        "effective_base_mass_exp_minus_divergence": math.exp(-divergence),
    }


def channel_arrays(
    left: np.ndarray, right: np.ndarray, total_order: int, orientation: int
) -> tuple[np.ndarray, np.ndarray, np.longdouble, dict]:
    rows, columns = len(left), len(right)
    dimension = rows * columns
    length = 1 << dimension
    t = np.longdouble(BETA) / np.sqrt(np.longdouble(total_order))
    x = exact.projective_spins(rows).astype(np.int16)
    y = exact.projective_spins(columns).astype(np.int16)
    ex = exact.energies_for_matrix(left, x)
    ey = exact.energies_for_matrix(right, y)

    weights = np.zeros(length, dtype=np.longdouble)
    patterns = []
    for xi, exi in zip(x, ex):
        for yj, eyj in zip(y, ey):
            signs = (xi[:, None] * yj[None, :]).reshape(-1)
            pattern = 0
            for bit, sign in enumerate(signs):
                if sign < 0:
                    pattern |= 1 << bit
            patterns.append(pattern)
            weights[pattern] += np.cosh(
                t * np.longdouble(exi + orientation * eyj)
            )
    if len(patterns) != len(set(patterns)):
        raise AssertionError("rank-one projective characters collided")

    indices = np.arange(length, dtype=np.uint64)
    popcount = np.bitwise_count(indices).astype(np.int16)
    overlap = dimension - 2 * popcount
    radial_t = np.cosh(t * overlap.astype(np.longdouble))
    radial_2t = np.cosh(2 * t * overlap.astype(np.longdouble))
    normalizer = np.longdouble(len(patterns))
    z_t = exact.xor_convolution(weights, radial_t) / normalizer
    z_2t = exact.xor_convolution(weights, radial_2t) / normalizer
    z_0 = np.sum(weights, dtype=np.longdouble) / normalizer
    if np.min(z_t) <= 0 or np.min(z_2t) <= 0 or z_0 <= 0:
        raise FloatingPointError("nonpositive channel pressure")

    # Direct checks use the same contracted-temperature internal weights but
    # bridge amplitudes t and 2t.  This guards against accidentally doubling
    # the child temperature along with the bridge temperature.
    check_masks = sorted(set((0, length // 3, length - 1)))
    maximum_log_errors = {"t": 0.0, "2t": 0.0}
    internal = ex[:, None] + orientation * ey[None, :]
    internal_weight = np.cosh(t * internal.astype(np.longdouble))
    for mask in check_masks:
        bits = ((mask >> np.arange(dimension, dtype=np.uint64)) & 1).astype(
            np.int16
        )
        bridge = (1 - 2 * bits).reshape(rows, columns)
        cross = x.astype(np.int64) @ bridge.astype(np.int64) @ y.astype(np.int64).T
        for label, amplitude, table in (
            ("t", t, z_t),
            ("2t", 2 * t, z_2t),
        ):
            direct = np.mean(
                internal_weight
                * np.cosh(amplitude * cross.astype(np.longdouble))
            )
            error = abs(float(np.log(direct) - np.log(table[mask])))
            maximum_log_errors[label] = max(maximum_log_errors[label], error)
    if max(maximum_log_errors.values()) > 3e-8:
        raise AssertionError(maximum_log_errors)
    return z_t, z_2t, z_0, {
        "bridge_dimension": dimension,
        "bridge_cube_size": length,
        "latent_projective_character_count": len(patterns),
        "direct_check_masks": check_masks,
        "maximum_direct_log_error": maximum_log_errors,
    }


def run_record(
    total_order: int,
    left: np.ndarray,
    right: np.ndarray,
    orientation: int,
) -> dict:
    started = time.time()
    z_t, z_2t, z_0, audit = channel_arrays(
        left, right, total_order, orientation
    )
    dimension = len(left) * len(right)
    length = len(z_t)
    t = BETA / math.sqrt(total_order)

    log_z_t = np.log(z_t).astype(np.float64)
    log_z_2t = np.log(z_2t).astype(np.float64)
    log_z_0 = float(np.log(z_0))
    shifted = -LAMBDA * log_z_t
    shifted -= float(np.max(shifted))
    escort = np.exp(shifted)
    escort /= float(np.sum(escort))

    log_k_full = log_z_0 + log_z_2t - 2.0 * log_z_t
    full_summary = collision_summary(log_k_full, escort, total_order)

    half = length // 2
    deleted_logs = np.empty(dimension * half, dtype=np.float64)
    deleted_weights = np.empty(dimension * half, dtype=np.float64)
    all_indices = np.arange(length, dtype=np.uint64)
    edge_means = []
    offset = 0
    for edge in range(dimension):
        bit = np.uint64(1 << edge)
        base = all_indices[(all_indices & bit) == 0]
        flipped = base | bit
        z_t_e = (z_t[base] + z_t[flipped]) / (
            2 * np.cosh(np.longdouble(t))
        )
        z_2t_e = (z_2t[base] + z_2t[flipped]) / (
            2 * np.cosh(np.longdouble(2 * t))
        )
        log_k_e = (
            log_z_0
            + np.log(z_2t_e).astype(np.float64)
            - 2.0 * np.log(z_t_e).astype(np.float64)
        )
        marginal_q = escort[base] + escort[flipped]
        deleted_logs[offset : offset + half] = log_k_e
        deleted_weights[offset : offset + half] = marginal_q / dimension
        edge_means.append(logsumexp_weighted(log_k_e, marginal_q))
        offset += half
    deleted_summary = collision_summary(
        deleted_logs, deleted_weights, total_order
    )
    deleted_summary["edge_log_mean_collision_over_N_range"] = [
        min(edge_means) / total_order,
        max(edge_means) / total_order,
    ]

    # The pointwise full/deleted comparison provides a useful numerical
    # regression check after independently computing both means.
    comparison_slack = (
        deleted_summary["log_mean_collision"]
        + 4 * t
        - full_summary["log_mean_collision"]
    )
    if comparison_slack < -2e-9:
        raise AssertionError("full/deleted comparison failed")

    return {
        "N": total_order,
        "split": [len(left), len(right)],
        "orientation": orientation,
        "beta": BETA,
        "lambda": LAMBDA,
        "raw_t": t,
        "channel_audit": audit,
        "inverse_escort": {
            "KL_from_fair_nats": float(
                np.dot(escort, np.log(escort * length))
            ),
            "entropy_nats": float(
                -np.dot(escort[escort > 0], np.log(escort[escort > 0]))
            ),
            "effective_bridge_words_exp_entropy": math.exp(
                float(-np.dot(escort[escort > 0], np.log(escort[escort > 0])))
            ),
        },
        "full_posterior_collision": full_summary,
        "deleted_posterior_collision": deleted_summary,
        "full_vs_deleted_log_comparison_slack": comparison_slack,
        "elapsed_seconds": time.time() - started,
    }


def main() -> None:
    mp.mp.dps = 80
    child_cache = {}
    records = []
    children = {}
    for total_order, left_order, right_order in PLANS:
        left, left_record = low_degree.child_record(
            left_order, total_order, child_cache
        )
        right, right_record = low_degree.child_record(
            right_order, total_order, child_cache
        )
        children[f"N{total_order}_left_{left_order}"] = left_record
        children[f"N{total_order}_right_{right_order}"] = right_record
        for orientation in ORIENTATIONS:
            record = run_record(total_order, left, right, orientation)
            records.append(record)
            full = record["full_posterior_collision"]
            deleted = record["deleted_posterior_collision"]
            print(
                f"N={total_order} split={left_order}+{right_order} "
                f"eps={orientation:+d} logK/N full="
                f"{full['log_mean_collision_over_N']:.9f} del="
                f"{deleted['log_mean_collision_over_N']:.9f} "
                f"elapsed={record['elapsed_seconds']:.1f}s",
                flush=True,
            )

    result = {
        "schema": "scratch-actual-child-latent-collision-v1",
        "classification": (
            "complete child/signing and bridge enumeration; interval child "
            "selection inherited from the certified routines; finite "
            "floating transcendental doubled-channel evaluation"
        ),
        "parameters": {
            "beta": BETA,
            "lambda": LAMBDA,
            "plans": [list(plan) for plan in PLANS],
            "orientations": list(ORIENTATIONS),
        },
        "identity": (
            "For fixed contracted-temperature child prior, K_full(B)="
            "z_2t(B) z_0 / z_t(B)^2; deleted K_e uses the corresponding "
            "edge-deleted z arrays. q_lambda is proportional to z_t^-lambda."
        ),
        "children": children,
        "records": records,
        "scope": (
            "finite complete-cube numerical evidence only; no asymptotic "
            "claim from N=6,8,10"
        ),
    }
    output = (
        ROOT
        / "computations/results/actual_child_latent_collision_diag.json"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
