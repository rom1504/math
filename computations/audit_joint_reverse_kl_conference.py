#!/usr/bin/env python3
"""Audit joint reverse-KL compensation on symmetric conference children.

For two copies of an order-r signing A, this program estimates (or, at
r=2, exhaustively evaluates)

    margin = 2 log Zbar_r(A, beta/sqrt(r))
             - E_{epsilon,B} log Zbar_{2r}([[A,B],[B^T,epsilon A]],
                                           beta/sqrt(2r)).

The proposed compensation inequality has the required direction only when
this margin is at least a geometrically summable negative error.  Child and
parent spin sums are exact for every sampled bridge.  The bridge average is
exhaustive at r=2 and seeded Monte Carlo at the larger orders.

No temporary files are used.  The saved exact signings in computations/results
are checked to satisfy A^2=(r-1)I before use.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def spins(order: int) -> np.ndarray:
    masks = np.arange(1 << order, dtype=np.uint64)[:, None]
    bits = (masks >> np.arange(order, dtype=np.uint64)) & 1
    return (1 - 2 * bits).astype(np.int16)


def conference(order: int) -> np.ndarray:
    if order == 2:
        matrix = np.asarray([[0, 1], [1, 0]], dtype=np.int16)
    elif order == 10:
        path = (
            ROOT
            / "computations"
            / "results"
            / "conference_order10_gf9.json"
        )
        matrix = np.asarray(
            json.loads(path.read_text())["conference_matrix"], dtype=np.int16
        )
    else:
        path = ROOT / "computations" / "results" / f"exact_m{order}.json"
        matrix = np.asarray(json.loads(path.read_text())["matrix"], dtype=np.int16)
    target = (order - 1) * np.eye(order, dtype=np.int16)
    if not np.array_equal(matrix @ matrix, target):
        raise ValueError(f"saved order-{order} signing is not conference")
    return matrix


def energies(matrix: np.ndarray, states: np.ndarray) -> np.ndarray:
    return (
        np.einsum("bi,ij,bj->b", states, matrix, states, dtype=np.int64) // 2
    )


def log_mean_cosh(energy: np.ndarray, raw_t: float, axis=None) -> np.ndarray:
    magnitude = np.abs(raw_t * np.asarray(energy, dtype=np.float64))
    peak = np.max(magnitude, axis=axis, keepdims=True)
    scaled = 0.5 * (
        np.exp(magnitude - peak) + np.exp(-magnitude - peak)
    )
    answer = peak + np.log(np.mean(scaled, axis=axis, keepdims=True))
    return np.squeeze(answer, axis=axis)


def parent_log_pressure(
    bridge: np.ndarray,
    orientation: int,
    states: np.ndarray,
    child_energy: np.ndarray,
    raw_t: float,
) -> float:
    cross = states @ bridge @ states.T
    total = (
        child_energy[:, None]
        + orientation * child_energy[None, :]
        + cross
    )
    return float(log_mean_cosh(total.reshape(-1), raw_t))


def audit_order(
    order: int,
    betas: list[float],
    samples: int,
    seed: int,
) -> dict:
    matrix = conference(order)
    states = spins(order)
    child_energy = energies(matrix, states)
    rng = np.random.default_rng(seed)

    if order == 2:
        bridge_count = 1 << (order * order)
        bridge_masks = range(bridge_count)
        status = "exhaustive bridge and orientation average"
    else:
        bridge_count = samples
        bridge_masks = range(samples)
        status = "seeded Monte Carlo bridge average; both orientations exact"

    sums = np.zeros(len(betas), dtype=np.float64)
    sum_squares = np.zeros(len(betas), dtype=np.float64)
    count = 0
    powers = np.arange(order * order, dtype=np.uint64)

    for mask in bridge_masks:
        if order == 2:
            bits = ((np.uint64(mask) >> powers) & 1).astype(np.int16)
            bridge = (1 - 2 * bits).reshape(order, order)
        else:
            bridge = rng.choice(
                np.asarray([-1, 1], dtype=np.int16), size=(order, order)
            )
        cross = states @ bridge @ states.T
        bridge_values = np.zeros(len(betas), dtype=np.float64)
        for orientation in (-1, 1):
            total = (
                child_energy[:, None]
                + orientation * child_energy[None, :]
                + cross
            ).reshape(-1)
            for index, beta in enumerate(betas):
                value = float(
                    log_mean_cosh(total, beta / math.sqrt(2 * order))
                )
                bridge_values[index] += 0.5 * value
        sums += bridge_values
        sum_squares += bridge_values * bridge_values
        count += 1

    records = []
    for index, beta in enumerate(betas):
        raw_t = beta / math.sqrt(2 * order)
        child_target = float(
            log_mean_cosh(child_energy, beta / math.sqrt(order))
        )
        child_low = float(log_mean_cosh(child_energy, raw_t))
        bridge_annealed = order * order * math.log(math.cosh(raw_t))
        mean_parent_log = sums[index] / count
        empirical_variance = max(
            0.0, sum_squares[index] / count - mean_parent_log**2
        )
        standard_error = (
            0.0 if order == 2 else math.sqrt(empirical_variance / (count - 1))
        )
        reverse_kl = 2 * child_low + bridge_annealed - mean_parent_log
        margin = 2 * child_target - mean_parent_log
        records.append(
            {
                "beta": beta,
                "raw_parent_temperature": raw_t,
                "child_target_log_Zbar": child_target,
                "child_contracted_log_Zbar": child_low,
                "bridge_annealed_cost": bridge_annealed,
                "mean_parent_log_Zbar": mean_parent_log,
                "mean_parent_log_standard_error": standard_error,
                "reverse_KL": reverse_kl,
                "compensation_margin": margin,
                "compensation_margin_per_parent_vertex": margin / (2 * order),
            }
        )

    beta4_coefficient = (
        -9 * order * order + 25 * order - 15
    ) / (48 * order)
    reverse_beta8_coefficient = 3 * (order - 1) ** 2 / (128 * order**2)
    return {
        "child_order": order,
        "parent_order": 2 * order,
        "matrix_square_verified": True,
        "bridge_draws": bridge_count,
        "orientation_values_per_bridge": 2,
        "output_samples": 2 * count,
        "independent_statistical_units": count,
        "seed": None if order == 2 else seed,
        "status": status,
        "proved_small_beta_coefficients": {
            "compensation_margin_beta2": -0.25,
            "compensation_margin_beta4": beta4_coefficient,
            "reverse_KL_beta8": reverse_beta8_coefficient,
        },
        "records": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orders", type=int, nargs="+", default=[2, 6, 10])
    parser.add_argument(
        "--betas", type=float, nargs="+", default=[0.25, 0.5, 1.0, 2.0]
    )
    parser.add_argument("--samples-6", type=int, default=8192)
    parser.add_argument("--samples-10", type=int, default=256)
    parser.add_argument("--seed", type=int, default=20260814)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    sample_counts = {6: args.samples_6, 10: args.samples_10}
    payload = {
        "schema": "joint-reverse-kl-conference-audit-v1",
        "classification": (
            "exact child and parent spin sums; exhaustive output average at "
            "order 2; seeded Monte Carlo output average at orders 6 and 10"
        ),
        "normalization": "Zbar_n=2^{-n} sum_x cosh(t H_A(x))",
        "orders": [
            audit_order(
                order,
                args.betas,
                sample_counts.get(order, 0),
                args.seed + order,
            )
            for order in args.orders
        ],
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered)


if __name__ == "__main__":
    main()
