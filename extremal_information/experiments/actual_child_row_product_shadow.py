#!/usr/bin/env python3
"""Coordinate-Gibbs audit of the actual child row-product variational shadow.

This is the preregistered companion to ``actual_child_bridge_law_exact.py``.
Every energy tensor is obtained by complete bridge enumeration.  Coordinate
updates solve the exact one-row subproblem, but global optimality is not
claimed because the product objective is nonconvex.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import mpmath as mp
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(HERE))
import actual_child_bridge_law_exact as exact  # noqa: E402


def softmax(log_weight: np.ndarray) -> np.ndarray:
    shifted = log_weight - np.max(log_weight)
    weight = np.exp(shifted)
    return weight / np.sum(weight)


def pressure_tensor(pressure: np.ndarray, rows: int, row_size: int) -> np.ndarray:
    # Flat bridge masks place row 0 in the least-significant block.  NumPy's
    # final reshape axis is least significant, hence reverse the axes.
    raw = pressure.reshape((row_size,) * rows)
    return np.transpose(raw, tuple(reversed(range(rows))))


def contract_except(
    tensor: np.ndarray,
    probabilities: list[np.ndarray],
    retained_axis: int | None,
) -> np.ndarray | float:
    current = tensor
    axes = list(range(len(probabilities)))
    for original_axis in reversed(range(len(probabilities))):
        if original_axis == retained_axis:
            continue
        current_axis = axes.index(original_axis)
        current = np.tensordot(
            current, probabilities[original_axis], axes=(current_axis, 0)
        )
        axes.pop(current_axis)
    if retained_axis is None:
        return float(current)
    if axes != [retained_axis]:
        raise AssertionError(axes)
    return np.asarray(current, dtype=np.float64)


def product_objective(
    tensor: np.ndarray,
    probabilities: list[np.ndarray],
    lam: float,
) -> float:
    energy = float(contract_except(tensor, probabilities, None))
    divergence = sum(
        exact.divergence_from_uniform(probability)
        for probability in probabilities
    )
    return energy + divergence / lam


def escort_row_marginal_start(
    pressure: np.ndarray,
    rows: int,
    columns: int,
    lam: float,
) -> list[np.ndarray]:
    q = softmax(-lam * pressure)
    return [
        exact.marginal_from_coordinates(q, exact.row_coordinates(row, columns))
        for row in range(rows)
    ]


def coordinate_descent(
    tensor: np.ndarray,
    initial: list[np.ndarray],
    lam: float,
    tolerance: float,
    max_sweeps: int,
) -> dict:
    probabilities = [probability.copy() for probability in initial]
    trace = [product_objective(tensor, probabilities, lam)]
    maximum_update_l1 = math.inf
    for sweep in range(1, max_sweeps + 1):
        maximum_update_l1 = 0.0
        for row in range(len(probabilities)):
            effective = contract_except(tensor, probabilities, row)
            updated = softmax(-lam * effective)
            maximum_update_l1 = max(
                maximum_update_l1,
                float(np.sum(np.abs(updated - probabilities[row]))),
            )
            probabilities[row] = updated
        objective = product_objective(tensor, probabilities, lam)
        if objective > trace[-1] + 5e-10:
            raise AssertionError((trace[-1], objective, sweep))
        trace.append(objective)
        if trace[-2] - trace[-1] <= tolerance:
            break

    # Recompute simultaneous best responses at the final point.  This is a
    # fixed-point residual, distinct from the last sequential update size.
    best_response_l1 = 0.0
    for row in range(len(probabilities)):
        effective = contract_except(tensor, probabilities, row)
        response = softmax(-lam * effective)
        best_response_l1 = max(
            best_response_l1,
            float(np.sum(np.abs(response - probabilities[row]))),
        )
    row_d2 = [
        math.log(len(probability) * float(np.dot(probability, probability)))
        for probability in probabilities
    ]
    return {
        "objective": trace[-1],
        "sweeps": len(trace) - 1,
        "converged_by_objective_tolerance": len(trace) - 1 < max_sweeps,
        "last_complete_sweep_decrease": trace[-2] - trace[-1],
        "last_sequential_maximum_update_l1": maximum_update_l1,
        "simultaneous_best_response_maximum_l1": best_response_l1,
        "row_Renyi2": row_d2,
        "objective_trace": trace,
    }


def optimize_product_shadow(
    pressure: np.ndarray,
    rows: int,
    columns: int,
    lam: float,
    random_starts: int,
    seed: int,
    tolerance: float,
    max_sweeps: int,
) -> dict:
    row_size = 1 << columns
    tensor = pressure_tensor(pressure, rows, row_size)
    rng = np.random.default_rng(seed)
    initializations: list[tuple[str, list[np.ndarray]]] = [
        (
            "uniform",
            [np.full(row_size, 1.0 / row_size) for _ in range(rows)],
        ),
        (
            "escort-row-marginals",
            escort_row_marginal_start(pressure, rows, columns, lam),
        ),
    ]
    scales = (0.1, 0.25, 0.5, 1.0, 2.0)
    for start in range(random_starts):
        scale = scales[start % len(scales)]
        initializations.append(
            (
                f"random-softmax-{start}-scale-{scale:g}",
                [softmax(scale * rng.standard_normal(row_size)) for _ in range(rows)],
            )
        )

    runs = []
    for name, initial in initializations:
        result = coordinate_descent(
            tensor, initial, lam, tolerance, max_sweeps
        )
        result["initialization"] = name
        runs.append(result)
    best = min(runs, key=lambda result: result["objective"])

    uniform_mean = float(np.mean(pressure))
    exact_soft = exact.negative_moment_soft_pressure(pressure, lam)
    total_gain = uniform_mean - exact_soft
    product_gain = uniform_mean - best["objective"]
    reverse_projection_upper = lam * (best["objective"] - exact_soft)
    if product_gain < -5e-9 or reverse_projection_upper < -5e-9:
        raise AssertionError((uniform_mean, exact_soft, best["objective"]))
    return {
        "lambda": lam,
        "uniform_mean_pressure": uniform_mean,
        "exact_negative_moment_soft_pressure": exact_soft,
        "best_evaluated_product_objective": best["objective"],
        "rigorous_candidate_row_product_gain": max(0.0, product_gain),
        "exact_total_escort_gain": total_gain,
        "captured_gain_fraction": (
            0.0 if total_gain <= 1e-15 else max(0.0, product_gain) / total_gain
        ),
        "candidate_reverse_projection_upper_bound": max(
            0.0, reverse_projection_upper
        ),
        "best_initialization": best["initialization"],
        "best_run": best,
        "all_start_summaries": [
            {
                key: value
                for key, value in run.items()
                if key != "objective_trace"
            }
            for run in runs
        ],
        "distinct_terminal_objective_count_1e-9": len(
            {round(run["objective"], 9) for run in runs}
        ),
    }


def audit(args: argparse.Namespace) -> dict:
    mp.mp.dps = args.mp_dps
    started = time.time()
    beta_texts = [format(beta, ".12g") for beta in args.betas]
    spaces = {
        n: exact.build_signing_space(n, args.signing_batch_size)
        for n in range(2, args.max_total_n - 1)
    }
    minimizer_cache: dict[tuple[int, str, int], tuple[list[dict], dict]] = {}
    records = []
    for total_n in range(args.min_total_n, args.max_total_n + 1):
        splits = list(range(2, total_n // 2 + 1))
        if total_n >= 9:
            splits = [total_n // 2]
        for m in splits:
            n = total_n - m
            for beta_text in beta_texts:
                beta = float(beta_text)
                class_lists = []
                for child_order in (m, n):
                    key = (child_order, beta_text, total_n)
                    if key not in minimizer_cache:
                        minimizer_cache[key] = exact.thermal_minimizer_classes(
                            spaces[child_order], beta_text, total_n
                        )
                    class_lists.append(minimizer_cache[key][0])
                for left_class in class_lists[0]:
                    left = np.asarray(
                        left_class["representative_matrix"], dtype=np.int8
                    )
                    for right_class in class_lists[1]:
                        right = np.asarray(
                            right_class["representative_matrix"], dtype=np.int8
                        )
                        for epsilon in (-1, 1):
                            pressure, pressure_audit = exact.bridge_pressures(
                                left, right, beta, total_n, epsilon
                            )
                            laws = []
                            for lambda_index, lam in enumerate(args.lambdas):
                                laws.append(
                                    optimize_product_shadow(
                                        pressure,
                                        m,
                                        n,
                                        lam,
                                        args.random_starts,
                                        args.seed
                                        + 100000 * total_n
                                        + 10000 * m
                                        + 1000 * int(round(100 * beta))
                                        + 100 * (epsilon + 1)
                                        + lambda_index,
                                        args.tolerance,
                                        args.max_sweeps,
                                    )
                                )
                            records.append(
                                {
                                    "N": total_n,
                                    "split": [m, n],
                                    "beta": beta,
                                    "raw_t": beta / math.sqrt(total_n),
                                    "relative_child_orientation": epsilon,
                                    "left_child_class": left_class["class_id"],
                                    "right_child_class": right_class["class_id"],
                                    "bridge_pressure_audit": pressure_audit,
                                    "laws": laws,
                                }
                            )
                            print(
                                f"N={total_n} split={m}+{n} beta={beta:g} "
                                f"eps={epsilon:+d} shadows={len(laws)}",
                                flush=True,
                            )
    return {
        "schema": "actual-child-row-product-shadow-coordinate-audit-v1",
        "classification": (
            "complete finite bridge enumeration plus nonconvex numerical "
            "coordinate minimization; candidate objectives are globally feasible "
            "but global row-product optimality is not certified"
        ),
        "protocol": "extremal_information/experiments/actual_child_bridge_law_protocol.md",
        "scope": {
            "N_le_8": "all splits",
            "N_9": "balanced 4+5 split only",
            "betas": args.betas,
            "lambdas": args.lambdas,
            "random_softmax_starts": args.random_starts,
            "additional_starts": ["uniform", "escort-row-marginals"],
            "seed": args.seed,
            "tolerance": args.tolerance,
            "max_sweeps": args.max_sweeps,
        },
        "records": records,
        "wall_time_seconds": time.time() - started,
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
    parser.add_argument("--random-starts", type=int, default=16)
    parser.add_argument("--seed", type=int, default=20260818)
    parser.add_argument("--tolerance", type=float, default=1e-11)
    parser.add_argument("--max-sweeps", type=int, default=100)
    parser.add_argument("--mp-dps", type=int, default=80)
    parser.add_argument("--signing-batch-size", type=int, default=8192)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT
        / "computations"
        / "results"
        / "actual_child_row_product_shadow.json",
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
