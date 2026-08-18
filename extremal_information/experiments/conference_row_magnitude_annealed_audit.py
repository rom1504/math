#!/usr/bin/env python3
"""Exact annealed audit for row-magnitude conditioned conference bridges.

For a conference child A of order r, fix a sign direction v and condition
each independent bridge row R on a sign-invariant event depending only on
|<R,v>|.  The full bridge family then has probability p_event**r, i.e. the
desired exp(-Theta(r)) scale when p_event is fixed away from zero and one.

This program computes E_B Zbar exactly without enumerating bridge matrices.
For a symmetric row law,

    E_B cosh(t(U + x^T B y)) = cosh(t U) M_v(y,t)^r,

where M_v(y,t)=E_R exp(t<R,y>) under the conditioned row law.  All row and
spin sums are exhaustive.  The result is annealed evidence only; it does not
identify E_B log Zbar or prove a quenched pressure theorem.

No temporary files are used.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]


def spins(order: int) -> np.ndarray:
    masks = np.arange(1 << order, dtype=np.uint64)[:, None]
    bits = (masks >> np.arange(order, dtype=np.uint64)) & 1
    return (1 - 2 * bits).astype(np.int16)


def conference(order: int) -> np.ndarray:
    if order == 10:
        payload = json.loads(
            (ROOT / "computations/results/conference_order10_gf9.json").read_text()
        )
        matrix = np.asarray(payload["conference_matrix"], dtype=np.int16)
    elif order == 14:
        payload = json.loads(
            (ROOT / "computations/results/conference_double_p13.json").read_text()
        )
        matrix = np.asarray(payload["conference_matrix"], dtype=np.int16)
    else:
        payload = json.loads(
            (ROOT / f"computations/results/exact_m{order}.json").read_text()
        )
        matrix = np.asarray(payload["matrix"], dtype=np.int16)
    target = (order - 1) * np.eye(order, dtype=np.int16)
    if not np.array_equal(matrix @ matrix, target):
        raise ValueError(f"order-{order} input is not conference")
    return matrix


def energies(matrix: np.ndarray, states: np.ndarray) -> np.ndarray:
    return np.einsum(
        "bi,ij,bj->b", states, matrix, states, dtype=np.int64
    ) // 2


def logsumexp(values: np.ndarray) -> float:
    peak = float(np.max(values))
    return peak + math.log(float(np.exp(values - peak).sum()))


def fwht(values: np.ndarray) -> np.ndarray:
    """Unnormalized in-place Walsh--Hadamard transform of a copy."""
    answer = np.asarray(values, dtype=np.float64).copy()
    width = 1
    size = answer.size
    while width < size:
        blocks = answer.reshape(-1, 2 * width)
        left = blocks[:, :width].copy()
        right = blocks[:, width:].copy()
        blocks[:, :width] = left + right
        blocks[:, width:] = left - right
        width *= 2
    return answer


def row_log_mgf(states: np.ndarray, row_mask: np.ndarray, t: float) -> np.ndarray:
    """All M(y)=E[exp(t<R,y>)|mask] by XOR convolution."""
    size, order = states.shape
    row_mass = row_mask.astype(np.float64) / float(row_mask.sum())
    kernel = np.exp(t * states.sum(axis=1).astype(np.float64))
    convolution = fwht(fwht(row_mass) * fwht(kernel)) / size
    if np.any(convolution <= 0):
        raise AssertionError("Walsh convolution produced a nonpositive MGF")
    return np.log(convolution)


def direction_bank(matrix: np.ndarray, states: np.ndarray, random_count: int, seed: int):
    order = matrix.shape[0]
    candidates: list[tuple[str, np.ndarray]] = [("all", np.ones(order, dtype=np.int16))]
    base = matrix + np.eye(order, dtype=np.int16)
    for index, row in enumerate(base):
        candidates.append((f"universal_double_row_{index}", row.astype(np.int16)))
    rng = np.random.default_rng(seed)
    for index in range(random_count):
        candidates.append(
            (f"random_{index}", rng.choice(np.asarray([-1, 1], dtype=np.int16), order))
        )
    if order <= 6:
        candidates.extend(
            (f"projective_{index}", state.copy())
            for index, state in enumerate(states[: 1 << (order - 1)])
        )
    unique: dict[bytes, tuple[str, np.ndarray]] = {}
    for label, vector in candidates:
        canonical = vector if vector[0] == 1 else -vector
        unique.setdefault(canonical.tobytes(), (label, canonical.copy()))
    return list(unique.values())


def event_masks(row_scores: np.ndarray) -> list[tuple[str, np.ndarray]]:
    magnitudes = np.abs(row_scores)
    levels = sorted(int(value) for value in np.unique(magnitudes))
    events: list[tuple[str, np.ndarray]] = []
    for level in levels[:-1]:
        events.append((f"abs_le_{level}", magnitudes <= level))
    for level in levels[1:]:
        events.append((f"abs_ge_{level}", magnitudes >= level))
    return events


def annealed_log_pressure(
    states: np.ndarray,
    energy: np.ndarray,
    row_mask: np.ndarray,
    beta: float,
    orientation: int,
) -> float:
    order = states.shape[1]
    t = beta / math.sqrt(2 * order)
    log_mgf = row_log_mgf(states, row_mask, t)

    scaled_energy = t * energy.astype(np.float64)
    c_child = float(np.mean(np.cosh(scaled_energy)))
    s_child = float(np.mean(np.sinh(scaled_energy)))
    y_factor = (
        c_child * np.cosh(scaled_energy)
        + orientation * s_child * np.sinh(scaled_energy)
    )
    if np.any(y_factor <= 0):
        raise AssertionError("positive conditional x-average became nonpositive")
    log_terms = order * log_mgf + np.log(y_factor)
    return logsumexp(log_terms) - math.log(states.shape[0])


def uniform_annealed_log_pressure(
    energy: np.ndarray, beta: float, orientation: int
) -> float:
    order = int(round(math.log2(energy.size)))
    t = beta / math.sqrt(2 * order)
    scaled = t * energy.astype(np.float64)
    c_child = float(np.mean(np.cosh(scaled)))
    s_child = float(np.mean(np.sinh(scaled)))
    internal = c_child * c_child + orientation * s_child * s_child
    return math.log(internal) + order * order * math.log(math.cosh(t))


def audit_order(order: int, betas: list[float], random_count: int, seed: int) -> dict:
    matrix = conference(order)
    states = spins(order)
    energy = energies(matrix, states)
    directions = direction_bank(matrix, states, random_count, seed + order)
    records = []
    for label, direction in directions:
        row_scores = states.astype(np.int64) @ direction.astype(np.int64)
        for event_label, mask in event_masks(row_scores):
            probability = float(mask.mean())
            if probability <= 0.0 or probability >= 1.0:
                continue
            for beta in betas:
                for orientation in (-1, 1):
                    value = annealed_log_pressure(
                        states, energy, mask, beta, orientation
                    )
                    baseline = uniform_annealed_log_pressure(
                        energy, beta, orientation
                    )
                    records.append(
                        {
                            "direction": label,
                            "event": event_label,
                            "row_event_probability": probability,
                            "full_bridge_log_probability_per_r": math.log(probability),
                            "beta": beta,
                            "orientation": orientation,
                            "conditioned_annealed_log_Zbar": value,
                            "uniform_annealed_log_Zbar": baseline,
                            "conditioned_minus_uniform": value - baseline,
                            "difference_per_r": (value - baseline) / order,
                        }
                    )
    records.sort(key=lambda item: item["difference_per_r"])
    nontrivial = [
        item for item in records
        if 0.1 <= item["row_event_probability"] <= 0.9
    ]
    minima_by_beta = {
        str(beta): min(
            item["difference_per_r"]
            for item in records
            if item["beta"] == beta
            and 0.1 <= item["row_event_probability"] <= 0.9
        )
        for beta in sorted({item["beta"] for item in records})
    }
    return {
        "order": order,
        "conference_sha256": hashlib.sha256(matrix.astype(np.int8).tobytes()).hexdigest(),
        "direction_count": len(directions),
        "directions": {
            label: vector.astype(int).tolist() for label, vector in directions
        },
        "record_count": len(records),
        "count_below_numerical_tolerance_minus_1e-8": sum(
            item["difference_per_r"] < -1e-8 for item in records
        ),
        "minimum_difference_per_r": records[0]["difference_per_r"],
        "minimum_nontrivial_probability_difference_per_r": nontrivial[0][
            "difference_per_r"
        ],
        "minimum_nontrivial_difference_by_beta": minima_by_beta,
        "smallest_conditioned_minus_uniform": records[: min(20, len(records))],
        "largest_conditioned_minus_uniform": records[-min(20, len(records)) :],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orders", nargs="+", type=int, default=[6, 10, 14])
    parser.add_argument("--betas", nargs="+", type=float, default=[0.1, 0.2, 0.5])
    parser.add_argument("--random-directions", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260818)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    payload = {
        "schema": "conference-row-magnitude-annealed-audit-v1",
        "classification": (
            "exact row and spin enumeration; exact annealed bridge pressure; "
            "not a quenched pressure computation or asymptotic theorem"
        ),
        "normalization": "Zbar=2^(-2r) sum_(x,y) cosh(beta H/sqrt(2r))",
        "orders": [
            audit_order(order, args.betas, args.random_directions, args.seed)
            for order in args.orders
        ],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(json.dumps({
        "output": str(args.output),
        "orders": args.orders,
        "record_counts": [item["record_count"] for item in payload["orders"]],
        "min_differences_per_r": [
            item["minimum_difference_per_r"] for item in payload["orders"]
        ],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
