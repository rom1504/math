#!/usr/bin/env python3
"""Exact finite audit of negative-disorder laws from actual thermal children.

The signing and bridge cubes are enumerated completely.  Integer energies and
energy histograms are exact.  Pressures, Gibbs weights, and information
quantities are numerical evaluations; child histogram comparisons are repeated
with mpmath at ``--mp-dps`` digits.

The preregistered methodology is in ``actual_child_bridge_law_protocol.md``.
No conference/Paley or ground-state surrogate is loaded by this program.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import mpmath as mp
import numpy as np


ROOT = Path(__file__).resolve().parents[2]
LN2 = math.log(2.0)


def projective_spins(n: int) -> np.ndarray:
    """All Boolean spins with the first spin fixed to +1."""
    if n == 0:
        return np.zeros((1, 0), dtype=np.int8)
    tails = np.asarray(
        list(itertools.product((-1, 1), repeat=max(0, n - 1))),
        dtype=np.int8,
    )
    return np.concatenate(
        [np.ones((len(tails), 1), dtype=np.int8), tails], axis=1
    )


def matrix_from_mask(n: int, mask: int) -> np.ndarray:
    matrix = np.zeros((n, n), dtype=np.int8)
    if n > 1:
        matrix[0, 1:] = 1
        matrix[1:, 0] = 1
    edges = [(i, j) for i in range(1, n) for j in range(i + 1, n)]
    for bit, (i, j) in enumerate(edges):
        value = -1 if (mask >> bit) & 1 else 1
        matrix[i, j] = matrix[j, i] = value
    return matrix


def rooted_gauge_key(matrix: np.ndarray) -> bytes:
    switches = np.ones(len(matrix), dtype=np.int8)
    if len(matrix) > 1:
        switches[1:] = matrix[0, 1:]
    gauged = switches[:, None] * matrix * switches[None, :]
    return gauged.astype(np.int8).tobytes()


def signed_permutation_orbit(reference: np.ndarray) -> set[bytes]:
    orbit: set[bytes] = set()
    n = len(reference)
    for permutation in itertools.permutations(range(n)):
        permuted = reference[np.ix_(permutation, permutation)]
        orbit.add(rooted_gauge_key(permuted))
        orbit.add(rooted_gauge_key(-permuted))
    return orbit


def matrix_sha(matrix: np.ndarray) -> str:
    return hashlib.sha256(matrix.astype(np.int8).tobytes()).hexdigest()


@dataclass
class SigningSpace:
    n: int
    spins: np.ndarray
    energies: np.ndarray
    absolute_histograms: np.ndarray
    histogram_inverse: np.ndarray
    unique_histograms: np.ndarray


def build_signing_space(n: int, batch_size: int = 8192) -> SigningSpace:
    spins = projective_spins(n).astype(np.int16)
    edges = [(i, j) for i in range(1, n) for j in range(i + 1, n)]
    total = 1 << len(edges)
    edge_products = np.asarray(
        [spins[:, i] * spins[:, j] for i, j in edges], dtype=np.int16
    ).T
    root_energy = np.sum(spins[:, 1:], axis=1, dtype=np.int16)
    energies = np.empty((total, len(spins)), dtype=np.int16)
    positions = np.arange(len(edges), dtype=np.uint64)
    for start in range(0, total, batch_size):
        masks = np.arange(start, min(total, start + batch_size), dtype=np.uint64)
        if edges:
            signs = 1 - 2 * (
                ((masks[:, None] >> positions) & 1).astype(np.int16)
            )
            energies[start : start + len(masks)] = (
                signs @ edge_products.T + root_energy[None, :]
            )
        else:
            energies[start : start + len(masks)] = root_energy[None, :]
    edge_count = n * (n - 1) // 2
    absolute = np.abs(energies)
    hist = np.stack(
        [np.count_nonzero(absolute == value, axis=1) for value in range(edge_count + 1)],
        axis=1,
    ).astype(np.int16)
    unique, inverse = np.unique(hist, axis=0, return_inverse=True)
    return SigningSpace(n, spins, energies, hist, inverse, unique)


def mp_hist_pressure(histogram: np.ndarray, beta_text: str, total_n: int) -> mp.mpf:
    t = mp.mpf(beta_text) / mp.sqrt(total_n)
    total = mp.fsum(
        mp.mpf(int(count)) * mp.cosh(t * value)
        for value, count in enumerate(histogram)
        if count
    )
    return mp.log(total / int(np.sum(histogram)))


def thermal_minimizer_classes(
    space: SigningSpace,
    beta_text: str,
    total_n: int,
) -> tuple[list[dict], dict]:
    """Return every signed-permutation class minimizing at beta/sqrt(total_n)."""
    t = float(beta_text) / math.sqrt(total_n)
    cosh_values = np.cosh(t * np.arange(space.unique_histograms.shape[1]))
    values = np.log(
        (space.unique_histograms.astype(np.float64) @ cosh_values)
        / len(space.spins)
    )
    order = np.argsort(values, kind="stable")
    # High precision is used on all histogram types.  At n<=7 their number is
    # small, and this avoids a hidden small-beta near-tie assumption.
    mp_values = [
        mp_hist_pressure(hist, beta_text, total_n)
        for hist in space.unique_histograms
    ]
    mp_optimum = min(mp_values)
    tie_tolerance = mp.mpf(10) ** (-(mp.mp.dps - 20))
    winner_hist_ids = [
        index
        for index, value in enumerate(mp_values)
        if abs(value - mp_optimum) <= tie_tolerance
    ]
    next_values = [
        value
        for index, value in enumerate(mp_values)
        if index not in winner_hist_ids
    ]
    next_value = min(next_values) if next_values else mp.inf
    winner_mask = np.isin(space.histogram_inverse, winner_hist_ids)
    masks = np.flatnonzero(winner_mask).tolist()

    key_to_mask = {
        rooted_gauge_key(matrix_from_mask(space.n, int(mask))): int(mask)
        for mask in masks
    }
    remaining = set(key_to_mask)
    classes: list[dict] = []
    while remaining:
        key = min(remaining)
        representative_mask = key_to_mask[key]
        representative = matrix_from_mask(space.n, representative_mask)
        orbit = signed_permutation_orbit(representative)
        members = sorted(key_to_mask[item] for item in remaining & orbit)
        if not members:
            raise AssertionError("orbit classification lost its representative")
        remaining.difference_update(orbit)
        classes.append(
            {
                "class_id": len(classes),
                "representative_mask": representative_mask,
                "representative_matrix": representative.astype(int).tolist(),
                "representative_sha256": matrix_sha(representative),
                "root_gauged_member_count": len(members),
                "root_gauged_member_masks": members,
            }
        )
    classes.sort(key=lambda row: row["representative_sha256"])
    for index, row in enumerate(classes):
        row["class_id"] = index

    double_gap = float(values[order[1]] - values[order[0]]) if len(order) > 1 else math.inf
    certificate = {
        "order": space.n,
        "total_parent_order": total_n,
        "beta": float(beta_text),
        "raw_t": t,
        "root_gauged_signing_count": int(len(space.energies)),
        "distinct_absolute_energy_histogram_count": int(len(space.unique_histograms)),
        "minimizing_histogram_count": len(winner_hist_ids),
        "minimizing_root_gauged_signing_count": len(masks),
        "signed_permutation_global_sign_class_count": len(classes),
        "mp_dps": mp.mp.dps,
        "mp_optimum": mp.nstr(mp_optimum, mp.mp.dps),
        "mp_gap_to_next_histogram": (
            None if next_value == mp.inf else mp.nstr(next_value - mp_optimum, mp.mp.dps)
        ),
        "double_nearest_sorted_histogram_gap": double_gap,
        "classification": (
            "exact signing/energy/histogram enumeration; high-precision numerical "
            "comparison of transcendental histogram pressures"
        ),
    }
    return classes, certificate


def energies_for_matrix(matrix: np.ndarray, spins: np.ndarray) -> np.ndarray:
    return (
        np.einsum(
            "bi,ij,bj->b",
            spins.astype(np.int64),
            matrix.astype(np.int64),
            spins.astype(np.int64),
            dtype=np.int64,
        )
        // 2
    )


def fwht(values: np.ndarray) -> np.ndarray:
    result = np.asarray(values, dtype=np.longdouble).copy()
    width = 1
    while width < len(result):
        view = result.reshape(-1, 2 * width)
        left = view[:, :width].copy()
        right = view[:, width:].copy()
        view[:, :width] = left + right
        view[:, width:] = left - right
        width *= 2
    return result


def xor_convolution(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    transformed = fwht(left) * fwht(right)
    return fwht(transformed) / np.longdouble(len(left))


def bridge_pressures(
    left: np.ndarray,
    right: np.ndarray,
    beta: float,
    total_n: int,
    epsilon: int,
) -> tuple[np.ndarray, dict]:
    m, n = len(left), len(right)
    d = m * n
    length = 1 << d
    t = np.longdouble(beta) / np.sqrt(np.longdouble(total_n))
    x = projective_spins(m).astype(np.int16)
    y = projective_spins(n).astype(np.int16)
    ex = energies_for_matrix(left, x)
    ey = energies_for_matrix(right, y)

    weights = np.zeros(length, dtype=np.longdouble)
    patterns: list[int] = []
    for xi, exi in zip(x, ex):
        for yj, eyj in zip(y, ey):
            signs = (xi[:, None] * yj[None, :]).reshape(-1)
            pattern = 0
            for bit, sign in enumerate(signs):
                if sign < 0:
                    pattern |= 1 << bit
            patterns.append(pattern)
            weights[pattern] += np.cosh(t * np.longdouble(exi + epsilon * eyj))
    if len(set(patterns)) != len(patterns):
        raise AssertionError("projective rank-one bridge characters collided")

    indices = np.arange(length, dtype=np.uint64)
    popcount = np.bitwise_count(indices).astype(np.uint8)
    radial = np.cosh(t * (d - 2 * popcount.astype(np.int16)))
    zbar = xor_convolution(weights, radial) / np.longdouble(len(patterns))
    minimum_z = float(np.min(zbar))
    if minimum_z <= 0:
        raise FloatingPointError(f"nonpositive partition value {minimum_z}")
    pressure = np.log(zbar).astype(np.float64)

    # Three deterministic direct checks guard the convolution/index convention.
    check_masks = sorted(set([0, length // 3, length - 1]))
    direct_errors: list[float] = []
    for mask in check_masks:
        bits = ((mask >> np.arange(d, dtype=np.uint64)) & 1).astype(np.int16)
        bridge = (1 - 2 * bits).reshape(m, n)
        cross = x.astype(np.int64) @ bridge.astype(np.int64) @ y.astype(np.int64).T
        internal = ex[:, None] + epsilon * ey[None, :]
        direct = np.mean(
            np.cosh(t * internal.astype(np.longdouble))
            * np.cosh(t * cross.astype(np.longdouble))
        )
        direct_errors.append(float(abs(np.log(direct) - pressure[mask])))
    max_error = max(direct_errors)
    if max_error > 2e-9:
        raise AssertionError((m, n, beta, epsilon, max_error))
    return pressure, {
        "bridge_sign_count": d,
        "bridge_cube_size": length,
        "rank_one_character_count": len(patterns),
        "convolution_dtype": str(np.dtype(np.longdouble)),
        "direct_check_masks": check_masks,
        "maximum_direct_log_pressure_error": max_error,
    }


def entropy(probability: np.ndarray) -> float:
    positive = probability[probability > 0]
    return float(-np.dot(positive, np.log(positive)))


def divergence_from_uniform(probability: np.ndarray) -> float:
    """Numerically stable D(p||uniform) near the uniform law."""
    positive = probability[probability > 0]
    return float(np.dot(positive, np.log(positive * len(probability))))


def weighted_quantile(values: np.ndarray, weights: np.ndarray, quantile: float) -> float:
    if len(values) == 0:
        return 0.0
    order = np.argsort(values, kind="stable")
    sorted_values = values[order]
    cumulative = np.cumsum(weights[order])
    target = quantile * cumulative[-1]
    index = min(len(values) - 1, int(np.searchsorted(cumulative, target, side="left")))
    return float(sorted_values[index])


def distribution_summary(
    values: np.ndarray,
    weights: np.ndarray,
    thresholds: Iterable[float] = (0.25, 0.5, 1.0, 2.0, 4.0, 8.0),
) -> dict:
    weights = np.asarray(weights, dtype=np.float64)
    weights = weights / np.sum(weights)
    values = np.asarray(values, dtype=np.float64)
    return {
        "weighted_mean": float(np.dot(weights, values)),
        "weighted_q50": weighted_quantile(values, weights, 0.50),
        "weighted_q90": weighted_quantile(values, weights, 0.90),
        "weighted_q99": weighted_quantile(values, weights, 0.99),
        "maximum": float(np.max(values)) if len(values) else 0.0,
        "mass_above_threshold": {
            format(threshold, ".12g"): float(np.sum(weights[values > threshold]))
            for threshold in thresholds
        },
    }


def marginal_from_coordinates(
    q: np.ndarray,
    coordinates: tuple[int, ...],
    cache: dict[tuple[int, ...], np.ndarray] | None = None,
) -> np.ndarray:
    if cache is not None and coordinates in cache:
        return cache[coordinates]
    if not coordinates:
        result = np.asarray([1.0], dtype=np.float64)
        if cache is not None:
            cache[coordinates] = result
        return result
    indices = np.arange(len(q), dtype=np.uint64)
    code = np.zeros(len(q), dtype=np.uint32)
    for target, source in enumerate(coordinates):
        code |= (((indices >> source) & 1).astype(np.uint32)) << target
    result = np.bincount(code, weights=q, minlength=1 << len(coordinates))
    if cache is not None:
        cache[coordinates] = result
    return result


def row_coordinates(row: int, row_width: int) -> tuple[int, ...]:
    return tuple(row * row_width + column for column in range(row_width))


def column_coordinates(column: int, rows: int, columns: int) -> tuple[int, ...]:
    return tuple(row * columns + column for row in range(rows))


def joint_block_distribution(
    q: np.ndarray,
    blocks: tuple[int, ...],
    block_coordinates: list[tuple[int, ...]],
    cache: dict[tuple[int, ...], np.ndarray] | None = None,
) -> np.ndarray:
    coordinates = tuple(
        coordinate for block in blocks for coordinate in block_coordinates[block]
    )
    return marginal_from_coordinates(q, coordinates, cache)


def conditional_d2_against_uniform(
    q: np.ndarray,
    latent_blocks: tuple[int, ...],
    target_block: int,
    block_coordinates: list[tuple[int, ...]],
    cache: dict[tuple[int, ...], np.ndarray] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    # The D2 distribution is unchanged by a permutation of the latent labels.
    latent_blocks = tuple(sorted(latent_blocks))
    width = len(block_coordinates[target_block])
    latent_coordinates = tuple(
        coordinate
        for block in latent_blocks
        for coordinate in block_coordinates[block]
    )
    target_coordinates = block_coordinates[target_block]
    joint = marginal_from_coordinates(
        q, latent_coordinates + target_coordinates, cache
    )
    latent_size = 1 << len(latent_coordinates)
    target_size = 1 << width
    table = joint.reshape(target_size, latent_size).T
    latent_probability = np.sum(table, axis=1)
    positive = latent_probability > 0
    conditional = table[positive] / latent_probability[positive, None]
    d2 = np.log(target_size * np.sum(conditional * conditional, axis=1))
    return d2, latent_probability[positive]


def block_information_metrics(
    q: np.ndarray,
    rows: int,
    columns: int,
) -> dict:
    row_blocks = [row_coordinates(row, columns) for row in range(rows)]
    column_blocks = [
        column_coordinates(column, rows, columns) for column in range(columns)
    ]
    hq = entropy(q)
    row_marginal_cache: dict[tuple[int, ...], np.ndarray] = {}
    column_marginal_cache: dict[tuple[int, ...], np.ndarray] = {}

    row_marginals = [
        joint_block_distribution(q, (i,), row_blocks, row_marginal_cache)
        for i in range(rows)
    ]
    row_entropies = [entropy(p) for p in row_marginals]
    row_d2 = [
        float(math.log((1 << columns) * float(np.dot(p, p))))
        for p in row_marginals
    ]
    row_tc = float(sum(row_entropies) - hq)

    column_marginals = [
        joint_block_distribution(q, (j,), column_blocks, column_marginal_cache)
        for j in range(columns)
    ]
    column_entropies = [entropy(p) for p in column_marginals]
    column_d2 = [
        float(math.log((1 << rows) * float(np.dot(p, p))))
        for p in column_marginals
    ]
    column_tc = float(sum(column_entropies) - hq)

    subset_entropy: dict[tuple[int, ...], float] = {(): 0.0}
    for size in range(1, rows + 1):
        for subset in itertools.combinations(range(rows), size):
            subset_entropy[subset] = entropy(
                joint_block_distribution(
                    q, subset, row_blocks, row_marginal_cache
                )
            )

    latent_records: list[dict] = []
    latent_choices = [()] + [(row,) for row in range(rows)]
    for latent in latent_choices:
        remaining = tuple(row for row in range(rows) if row not in latent)
        hz = subset_entropy[tuple(sorted(latent))]
        if len(remaining) <= 1:
            residual_tc = 0.0
        else:
            conditional_entropy_sum = 0.0
            for row in remaining:
                union = tuple(sorted(latent + (row,)))
                conditional_entropy_sum += subset_entropy[union] - hz
            residual_tc = conditional_entropy_sum - (hq - hz)
        all_values: list[np.ndarray] = []
        all_weights: list[np.ndarray] = []
        per_target: list[dict] = []
        for row in remaining:
            d2_values, prefix_weights = conditional_d2_against_uniform(
                q, latent, row, row_blocks, row_marginal_cache
            )
            per_target.append(
                {
                    "target_row": row,
                    **distribution_summary(d2_values, prefix_weights),
                }
            )
            all_values.append(d2_values)
            all_weights.append(prefix_weights / max(1, len(remaining)))
        aggregate = (
            distribution_summary(np.concatenate(all_values), np.concatenate(all_weights))
            if all_values
            else distribution_summary(np.asarray([0.0]), np.asarray([1.0]))
        )
        latent_records.append(
            {
                "latent_rows": list(latent),
                "latent_state_log_cardinality": len(latent) * columns * LN2,
                "remaining_row_count": len(remaining),
                "residual_conditional_total_correlation": max(0.0, float(residual_tc)),
                "conditional_marginal_d2_aggregate": aggregate,
                "conditional_marginal_d2_by_target": per_target,
            }
        )

    filtration_records: list[dict] = []
    for permutation in itertools.permutations(range(rows)):
        step_records: list[dict] = []
        all_values = []
        all_weights = []
        for position, row in enumerate(permutation):
            prefix = tuple(permutation[:position])
            values, weights = conditional_d2_against_uniform(
                q, prefix, row, row_blocks, row_marginal_cache
            )
            step_records.append(
                {
                    "position": position,
                    "row": row,
                    "prefix_rows": list(prefix),
                    **distribution_summary(values, weights),
                }
            )
            all_values.append(values)
            all_weights.append(weights / rows)
        aggregate = distribution_summary(
            np.concatenate(all_values), np.concatenate(all_weights)
        )
        filtration_records.append(
            {
                "row_order": list(permutation),
                "aggregate": aggregate,
                "maximum_step_weighted_mean": max(
                    step["weighted_mean"] for step in step_records
                ),
                "maximum_step_mass_above_1": max(
                    step["mass_above_threshold"]["1"] for step in step_records
                ),
                "_steps": step_records,
            }
        )
    best_filtration = min(
        filtration_records,
        key=lambda row: (
            row["maximum_step_weighted_mean"],
            row["aggregate"]["weighted_mean"],
        ),
    )
    best_latent_tc = min(
        latent_records,
        key=lambda row: row["residual_conditional_total_correlation"],
    )
    best_latent_d2 = min(
        latent_records,
        key=lambda row: row["conditional_marginal_d2_aggregate"]["weighted_mean"],
    )
    return {
        "row_total_correlation": max(0.0, row_tc),
        "row_total_correlation_per_row": max(0.0, row_tc) / rows,
        "row_marginal_d2": row_d2,
        "column_total_correlation": max(0.0, column_tc),
        "column_total_correlation_per_column": max(0.0, column_tc) / columns,
        "column_marginal_d2": column_d2,
        "latent_product_audits": latent_records,
        "best_latent_residual_tc": {
            "latent_rows": best_latent_tc["latent_rows"],
            "value": best_latent_tc["residual_conditional_total_correlation"],
        },
        "best_latent_conditional_d2_mean": {
            "latent_rows": best_latent_d2["latent_rows"],
            "value": best_latent_d2["conditional_marginal_d2_aggregate"][
                "weighted_mean"
            ],
        },
        "best_row_filtration": {
            "row_order": best_filtration["row_order"],
            "maximum_step_weighted_mean": best_filtration[
                "maximum_step_weighted_mean"
            ],
            "aggregate": best_filtration["aggregate"],
            "steps": best_filtration["_steps"],
        },
        "row_filtrations": [
            {key: value for key, value in record.items() if key != "_steps"}
            for record in filtration_records
        ],
    }


def natural_coordinate_orders(rows: int, columns: int, q: np.ndarray) -> list[tuple[str, tuple[int, ...]]]:
    row = tuple(range(rows * columns))
    column = tuple(i * columns + j for j in range(columns) for i in range(rows))
    row_snake = tuple(
        i * columns + j
        for i in range(rows)
        for j in (range(columns) if i % 2 == 0 else reversed(range(columns)))
    )
    column_snake = tuple(
        i * columns + j
        for j in range(columns)
        for i in (range(rows) if j % 2 == 0 else reversed(range(rows)))
    )
    marginal_kl: list[tuple[float, int]] = []
    for coordinate in range(rows * columns):
        p = marginal_from_coordinates(q, (coordinate,))
        marginal_kl.append((divergence_from_uniform(p), coordinate))
    ascending = tuple(coordinate for _, coordinate in sorted(marginal_kl))
    descending = tuple(reversed(ascending))
    candidates = [
        ("row-major", row),
        ("row-major-reverse", tuple(reversed(row))),
        ("column-major", column),
        ("column-major-reverse", tuple(reversed(column))),
        ("row-snake", row_snake),
        ("column-snake", column_snake),
        ("marginal-KL-ascending", ascending),
        ("marginal-KL-descending", descending),
    ]
    unique: list[tuple[str, tuple[int, ...]]] = []
    seen: set[tuple[int, ...]] = set()
    for name, order in candidates:
        if order not in seen:
            unique.append((name, order))
            seen.add(order)
    return unique


def reindex_for_order(d: int, order: tuple[int, ...]) -> np.ndarray:
    ordered = np.arange(1 << d, dtype=np.uint64)
    natural = np.zeros(1 << d, dtype=np.uint64)
    for target, source in enumerate(order):
        natural |= ((ordered >> target) & 1) << source
    return natural.astype(np.int64)


def chain_support_for_order(q: np.ndarray, order: tuple[int, ...]) -> tuple[float, list[float]]:
    d = len(order)
    p = q[reindex_for_order(d, order)]
    prefix_d = np.zeros(d + 1, dtype=np.float64)
    prefix_d[d] = divergence_from_uniform(p)
    current = p
    for j in range(d - 1, 0, -1):
        current = current.reshape(2, -1).sum(axis=0)
        prefix_d[j] = divergence_from_uniform(current)
    total = prefix_d[d]
    increments = np.maximum(0.0, np.diff(prefix_d))
    increments[increments <= max(1e-14, abs(total) * 1e-11)] = 0.0
    support = 0.0 if total <= 1e-15 else float(np.sum(np.sqrt(increments)) ** 2 / total)
    return support, increments.tolist()


def exact_chain_support(q: np.ndarray, d: int) -> float | None:
    if d > 10:
        return None
    length = 1 << d
    indices = np.arange(length, dtype=np.uint32)
    divergence = np.zeros(length, dtype=np.float64)
    for subset in range(1, length):
        coordinates = [j for j in range(d) if (subset >> j) & 1]
        code = np.zeros(length, dtype=np.uint16)
        for target, source in enumerate(coordinates):
            code |= (((indices >> source) & 1).astype(np.uint16)) << target
        marginal = np.bincount(code, weights=q, minlength=1 << len(coordinates))
        divergence[subset] = divergence_from_uniform(marginal)
    dp = np.full(length, np.inf, dtype=np.float64)
    dp[0] = 0.0
    for subset in range(1, length):
        bits = subset
        while bits:
            bit = bits & -bits
            previous = subset ^ bit
            increment = max(0.0, divergence[subset] - divergence[previous])
            if increment <= max(1e-14, abs(divergence[-1]) * 1e-11):
                increment = 0.0
            dp[subset] = min(dp[subset], dp[previous] + math.sqrt(increment))
            bits ^= bit
    total = divergence[-1]
    return 0.0 if total <= 1e-15 else float(dp[-1] ** 2 / total)


def law_metrics(
    pressure: np.ndarray,
    lam: float,
    rows: int,
    columns: int,
    total_n: int,
    same_temperature_child_target: float,
) -> dict:
    log_weight = -lam * pressure
    peak = float(np.max(log_weight))
    unnormalized = np.exp(log_weight - peak)
    q = unnormalized / np.sum(unnormalized)
    log_mean_exponential = peak + math.log(float(np.mean(unnormalized)))
    positive_peak = float(np.max(pressure))
    positive_unnormalized = np.exp(pressure - positive_peak)
    p = positive_unnormalized / np.sum(positive_unnormalized)
    log_mean_positive = positive_peak + math.log(
        float(np.mean(positive_unnormalized))
    )
    mean_uniform = float(np.mean(pressure))
    mean_q = float(np.dot(q, pressure))
    mean_p = float(np.dot(p, pressure))
    soft_pressure = -log_mean_exponential / lam
    hq = entropy(q)
    d = rows * columns
    kl = divergence_from_uniform(q)
    kl_identity = -lam * mean_q - log_mean_exponential
    d2_full = math.log((1 << d) * float(np.dot(q, q)))
    d_uniform_p = log_mean_positive - mean_uniform
    d_q_p = (-lam - 1.0) * mean_q - log_mean_exponential + log_mean_positive
    d_p_q = (lam + 1.0) * mean_p + log_mean_exponential - log_mean_positive
    tv_q_p = 0.5 * float(np.sum(np.abs(q - p)))
    affinity_q_p = float(np.sum(np.sqrt(q * p)))
    chain_records = []
    for name, order in natural_coordinate_orders(rows, columns, q):
        support, increments = chain_support_for_order(q, order)
        chain_records.append(
            {
                "ordering": name,
                "coordinate_order": list(order),
                "effective_support_proxy": support,
                "conditional_kl_increments": increments,
            }
        )
    best_chain = min(chain_records, key=lambda row: row["effective_support_proxy"])
    exact_support_raw = exact_chain_support(q, d)
    exact_support_correction = 0.0
    if (
        exact_support_raw is not None
        and exact_support_raw > best_chain["effective_support_proxy"]
    ):
        # The subset DP and the displayed chain are mathematically ordered in
        # the opposite direction.  Near uniformity, subtracting marginal
        # entropies loses a few digits; clamp only this certified numerical
        # inconsistency and retain its size for audit.
        exact_support_correction = (
            exact_support_raw - best_chain["effective_support_proxy"]
        )
        if exact_support_correction > 1e-4:
            raise AssertionError(exact_support_correction)
        exact_support = best_chain["effective_support_proxy"]
    else:
        exact_support = exact_support_raw
    block = block_information_metrics(q, rows, columns)
    return {
        "lambda": lam,
        "uniform_mean_pressure": mean_uniform,
        "tilted_mean_pressure": mean_q,
        "mean_pressure_gain": mean_uniform - mean_q,
        "mean_pressure_gain_per_parent_vertex": (mean_uniform - mean_q) / total_n,
        "negative_moment_soft_pressure": soft_pressure,
        "negative_moment_pressure_gain": mean_uniform + log_mean_exponential / lam,
        "negative_moment_pressure_gain_per_parent_vertex": (
            mean_uniform + log_mean_exponential / lam
        )
        / total_n,
        "same_temperature_child_target": same_temperature_child_target,
        "tilted_mean_gap_above_child_target": mean_q
        - same_temperature_child_target,
        "tilted_mean_gap_above_child_target_per_parent_vertex": (
            mean_q - same_temperature_child_target
        )
        / total_n,
        "negative_moment_gap_above_child_target": soft_pressure
        - same_temperature_child_target,
        "negative_moment_gap_above_child_target_per_parent_vertex": (
            soft_pressure - same_temperature_child_target
        )
        / total_n,
        "KL_q_parallel_U": kl,
        "KL_identity_residual": kl - kl_identity,
        "KL_per_parent_vertex": kl / total_n,
        "KL_per_bridge_sign": kl / d,
        "Renyi2_q_parallel_U": d2_full,
        "positive_output_law_mean_pressure": mean_p,
        "KL_U_parallel_positive_output_law": d_uniform_p,
        "KL_U_parallel_positive_output_law_per_parent_vertex": d_uniform_p
        / total_n,
        "KL_q_parallel_positive_output_law": d_q_p,
        "KL_q_parallel_positive_output_law_per_parent_vertex": d_q_p
        / total_n,
        "KL_positive_output_law_parallel_q": d_p_q,
        "KL_positive_output_law_parallel_q_per_parent_vertex": d_p_q
        / total_n,
        "TV_q_positive_output_law": tv_q_p,
        "affinity_q_positive_output_law": affinity_q_p,
        "Shannon_effective_bridge_count_log": hq,
        "Shannon_effective_bridge_fraction_log": -kl,
        "exact_chain_effective_support": exact_support,
        "exact_chain_effective_support_raw": exact_support_raw,
        "exact_chain_numerical_order_correction": exact_support_correction,
        "chain_support_proxy_minimum": best_chain["effective_support_proxy"],
        "chain_support_proxy_best_ordering": best_chain["ordering"],
        "chain_support_proxy_fraction": best_chain["effective_support_proxy"] / d,
        "chain_ordering_audits": chain_records,
        **block,
    }


def negative_moment_soft_pressure(pressure: np.ndarray, lam: float) -> float:
    log_weight = -lam * pressure
    peak = float(np.max(log_weight))
    return -(
        peak + math.log(float(np.mean(np.exp(log_weight - peak))))
    ) / lam


def target_threshold_lambda(
    pressure: np.ndarray,
    target: float,
    tolerance: float = 1e-11,
) -> dict:
    """Find the least disorder inverse temperature whose soft pressure hits target."""
    uniform_mean = float(np.mean(pressure))
    minimum = float(np.min(pressure))
    if minimum > target + tolerance:
        return {
            "status": "impossible-even-at-zero-disorder-temperature",
            "minimum_bridge_pressure_gap": minimum - target,
            "lambda": None,
        }
    if uniform_mean <= target + tolerance:
        return {
            "status": "uniform-law-already-at-target",
            "minimum_bridge_pressure_gap": minimum - target,
            "lambda": 0.0,
        }
    lower = 0.0
    upper = 1.0
    while negative_moment_soft_pressure(pressure, upper) > target:
        lower = upper
        upper *= 2.0
        if upper > 65536:
            return {
                "status": "numerical-threshold-not-bracketed",
                "minimum_bridge_pressure_gap": minimum - target,
                "lambda": None,
            }
    for _ in range(70):
        middle = 0.5 * (lower + upper)
        if negative_moment_soft_pressure(pressure, middle) > target:
            lower = middle
        else:
            upper = middle
    return {
        "status": "finite-target-threshold",
        "minimum_bridge_pressure_gap": minimum - target,
        "lambda": upper,
        "bracket_width": upper - lower,
        "soft_pressure_residual": negative_moment_soft_pressure(pressure, upper)
        - target,
    }


def audit(args: argparse.Namespace) -> dict:
    mp.mp.dps = args.mp_dps
    start_time = time.time()
    beta_texts = [format(value, ".12g") for value in args.betas]
    spaces = {
        n: build_signing_space(n, args.signing_batch_size)
        for n in range(2, args.max_total_n - 1)
    }
    minimizer_cache: dict[tuple[int, str, int], tuple[list[dict], dict]] = {}

    records: list[dict] = []
    child_certificates: list[dict] = []
    for total_n in range(args.min_total_n, args.max_total_n + 1):
        for m in range(2, total_n // 2 + 1):
            n = total_n - m
            for beta_text in beta_texts:
                beta = float(beta_text)
                class_lists = []
                for child_order in (m, n):
                    key = (child_order, beta_text, total_n)
                    if key not in minimizer_cache:
                        minimizer_cache[key] = thermal_minimizer_classes(
                            spaces[child_order], beta_text, total_n
                        )
                        child_certificates.append(minimizer_cache[key][1])
                    class_lists.append(minimizer_cache[key][0])
                left_classes, right_classes = class_lists
                for left_class in left_classes:
                    left = np.asarray(left_class["representative_matrix"], dtype=np.int8)
                    for right_class in right_classes:
                        right = np.asarray(
                            right_class["representative_matrix"], dtype=np.int8
                        )
                        for epsilon in (-1, 1):
                            pressure, pressure_audit = bridge_pressures(
                                left, right, beta, total_n, epsilon
                            )
                            left_contracted = float(
                                minimizer_cache[(m, beta_text, total_n)][1][
                                    "mp_optimum"
                                ]
                            )
                            right_contracted = float(
                                minimizer_cache[(n, beta_text, total_n)][1][
                                    "mp_optimum"
                                ]
                            )
                            target_certificates = []
                            for child_order in (m, n):
                                target_key = (child_order, beta_text, child_order)
                                if target_key not in minimizer_cache:
                                    minimizer_cache[target_key] = thermal_minimizer_classes(
                                        spaces[child_order], beta_text, child_order
                                    )
                                    child_certificates.append(
                                        minimizer_cache[target_key][1]
                                    )
                                target_certificates.append(
                                    minimizer_cache[target_key][1]
                                )
                            same_temperature_target = sum(
                                float(certificate["mp_optimum"])
                                for certificate in target_certificates
                            )
                            laws = [
                                law_metrics(
                                    pressure,
                                    lam,
                                    m,
                                    n,
                                    total_n,
                                    same_temperature_target,
                                )
                                for lam in args.lambdas
                            ]
                            threshold = target_threshold_lambda(
                                pressure, same_temperature_target
                            )
                            if threshold["status"] == "finite-target-threshold":
                                threshold["law"] = law_metrics(
                                    pressure,
                                    float(threshold["lambda"]),
                                    m,
                                    n,
                                    total_n,
                                    same_temperature_target,
                                )
                            records.append(
                                {
                                    "N": total_n,
                                    "split": [m, n],
                                    "beta": beta,
                                    "raw_t": beta / math.sqrt(total_n),
                                    "left_child_class": left_class["class_id"],
                                    "right_child_class": right_class["class_id"],
                                    "left_child_representative_sha256": left_class[
                                        "representative_sha256"
                                    ],
                                    "right_child_representative_sha256": right_class[
                                        "representative_sha256"
                                    ],
                                    "relative_child_orientation": epsilon,
                                    "contracted_temperature_child_pressure_sum": (
                                        left_contracted + right_contracted
                                    ),
                                    "same_temperature_minimum_child_target": (
                                        same_temperature_target
                                    ),
                                    "child_thermal_gap": (
                                        same_temperature_target
                                        - left_contracted
                                        - right_contracted
                                    ),
                                    "bridge_pressure_audit": pressure_audit,
                                    "bridge_pressure_minimum": float(np.min(pressure)),
                                    "bridge_pressure_maximum": float(np.max(pressure)),
                                    "bridge_pressure_standard_deviation": float(
                                        np.std(pressure)
                                    ),
                                    "target_threshold": threshold,
                                    "laws": laws,
                                }
                            )
                            print(
                                f"N={total_n} split={m}+{n} beta={beta:g} "
                                f"classes={left_class['class_id']},{right_class['class_id']} "
                                f"eps={epsilon:+d} cube=2^{m*n}",
                                flush=True,
                            )
    return {
        "schema": "actual-child-negative-disorder-bridge-law-audit-v1",
        "classification": (
            "exact finite signing and bridge enumeration; integer-exact energies; "
            "high-precision numerical thermal-child selection; floating-point "
            "Gibbs and information metrics"
        ),
        "protocol": "extremal_information/experiments/actual_child_bridge_law_protocol.md",
        "normalization": (
            "p_t(A)=log(2^-n sum_x cosh(t H_A(x))); "
            "q_lambda(B) proportional exp(-lambda f_A,C,epsilon(B))"
        ),
        "scope": {
            "min_total_order": args.min_total_n,
            "max_total_order": args.max_total_n,
            "splits": "all 2<=m<=n, m+n=N",
            "betas": args.betas,
            "lambdas": args.lambdas,
            "relative_child_orientations": [-1, 1],
            "mp_dps": args.mp_dps,
        },
        "exclusions": [
            "no conference/Paley child",
            "no ground-state child surrogate",
            "no bridge sampling",
            "no pressure-threshold-defined latent state",
        ],
        "child_minimizer_certificates": child_certificates,
        "records": records,
        "wall_time_seconds": time.time() - start_time,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-total-n", type=int, default=4)
    parser.add_argument("--max-total-n", type=int, default=9)
    parser.add_argument(
        "--betas", type=float, nargs="+", default=[0.25, 0.5, 1.0, 2.0, 4.0]
    )
    parser.add_argument(
        "--lambdas", type=float, nargs="+", default=[0.25, 0.5, 1.0, 2.0, 4.0]
    )
    parser.add_argument("--mp-dps", type=int, default=80)
    parser.add_argument("--signing-batch-size", type=int, default=8192)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT
        / "computations"
        / "results"
        / "actual_child_bridge_law_exact.json",
    )
    args = parser.parse_args()
    payload = audit(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
