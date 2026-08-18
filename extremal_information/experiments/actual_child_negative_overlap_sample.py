#!/usr/bin/env python3
"""Reproducible held-out sampling of actual-child negative overlap.

This complements ``actual_child_negative_overlap_exact.py`` beyond the exact
bridge-cube range.  Child minimizers are still selected by exhaustive signing
enumeration.  Only fair bridges are sampled.  For each sampled bridge the
full finite Gibbs sum is evaluated, and the exact full-to-cavity Möbius
identity converts every edge response to ``r_e(B_-e)``.  Self-normalized
importance sampling then evaluates the raw negative-tilt path.

The output is finite numerical evidence, not a certificate or an asymptotic
claim.  The default N=8 sample is compared with the complete bridge-cube
enumeration as an implementation check.
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


def sample_pressure_and_cavity(
    left_matrix: np.ndarray,
    right_matrix: np.ndarray,
    beta: float,
    total_order: int,
    orientation: int,
    sample_count: int,
    batch_size: int,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample ``L(B)`` and ``(mn)^-1 sum_e r_e(B_-e)^2``."""

    rows = len(left_matrix)
    columns = len(right_matrix)
    edge_count = rows * columns
    raw_t = beta / math.sqrt(total_order)
    channel = math.tanh(raw_t)
    x = exact.projective_spins(rows).astype(np.int16)
    y = exact.projective_spins(columns).astype(np.int16)
    ex = exact.energies_for_matrix(left_matrix, x)
    ey = exact.energies_for_matrix(right_matrix, y)
    patterns = []
    internal_weights = []
    for left_spin, left_energy in zip(x, ex):
        for right_spin, right_energy in zip(y, ey):
            patterns.append(
                (left_spin[:, None] * right_spin[None, :]).reshape(-1)
            )
            internal_weights.append(
                math.cosh(
                    raw_t
                    * float(left_energy + orientation * right_energy)
                )
            )
    rank_one = np.asarray(patterns, dtype=np.float64)
    internal = np.asarray(internal_weights, dtype=np.float64)

    pressures = []
    cavity_densities = []
    signs = np.asarray([-1.0, 1.0])
    for start in range(0, sample_count, batch_size):
        size = min(batch_size, sample_count - start)
        bridge = rng.choice(signs, size=(size, edge_count))
        field = bridge @ rank_one.T
        cosh_field = np.cosh(raw_t * field)
        sinh_field = np.sinh(raw_t * field)
        denominator = cosh_field @ internal
        full_response = (
            (sinh_field * internal[None, :]) @ rank_one
        ) / denominator[:, None]
        cavity_response = (
            full_response - channel * bridge
        ) / (1.0 - channel * bridge * full_response)
        if float(np.max(np.abs(cavity_response))) > 1.0 + 2e-10:
            raise AssertionError(float(np.max(np.abs(cavity_response))))
        pressures.append(np.log(denominator / len(rank_one)))
        cavity_densities.append(
            np.mean(cavity_response * cavity_response, axis=1)
        )
    return np.concatenate(pressures), np.concatenate(cavity_densities)


def ratio_with_influence(
    pressure: np.ndarray, values: np.ndarray, tilt: float
) -> tuple[float, np.ndarray, float]:
    exponent = tilt * pressure
    exponent -= float(np.max(exponent))
    weight = np.exp(exponent)
    mean_weight = float(np.mean(weight))
    mean = float(np.mean(weight * values) / mean_weight)
    influence = weight * (values - mean) / mean_weight
    effective_fraction = float(
        np.sum(weight) ** 2 / (len(weight) * np.dot(weight, weight))
    )
    return mean, influence, effective_fraction


def negative_path(
    pressure: np.ndarray,
    values: np.ndarray,
    inverse_disorder: float,
    quadrature_order: int,
) -> dict:
    nodes, quadrature_weights = np.polynomial.legendre.leggauss(
        quadrature_order
    )
    path_mean = 0.0
    path_influence = np.zeros(len(values), dtype=np.float64)
    minimum_effective_fraction = 1.0
    for node, quadrature_weight in zip(nodes, quadrature_weights):
        tilt = 0.5 * inverse_disorder * (float(node) - 1.0)
        mean, influence, effective_fraction = ratio_with_influence(
            pressure, values, tilt
        )
        coefficient = 0.5 * float(quadrature_weight)
        path_mean += coefficient * mean
        path_influence += coefficient * influence
        minimum_effective_fraction = min(
            minimum_effective_fraction, effective_fraction
        )
    endpoint, endpoint_influence, endpoint_effective_fraction = (
        ratio_with_influence(pressure, values, -inverse_disorder)
    )
    uniform, uniform_influence, _ = ratio_with_influence(
        pressure, values, 0.0
    )

    def standard_error(influence: np.ndarray) -> float:
        return float(np.std(influence, ddof=1) / math.sqrt(len(influence)))

    return {
        "lambda": inverse_disorder,
        "uniform_overlap": uniform,
        "uniform_overlap_standard_error": standard_error(uniform_influence),
        "rhohat_path_average": path_mean,
        "rhohat_path_standard_error": standard_error(path_influence),
        "endpoint_overlap_at_minus_lambda": endpoint,
        "endpoint_overlap_standard_error": standard_error(
            endpoint_influence
        ),
        "minimum_quadrature_importance_effective_fraction": (
            minimum_effective_fraction
        ),
        "endpoint_importance_effective_fraction": (
            endpoint_effective_fraction
        ),
    }


def exact_reference(
    path: Path,
    total_order: int,
    beta: float,
    orientation: int,
    inverse_disorder: float,
) -> dict | None:
    if not path.exists():
        return None
    payload = json.loads(path.read_text(encoding="utf-8"))
    for record in payload["records"]:
        if (
            record["N"] == total_order
            and record["beta"] == beta
            and record["orientation"] == orientation
            and record["left_child_class"] == 0
            and record["right_child_class"] == 0
        ):
            for candidate in record["paths"]:
                if candidate["lambda"] == inverse_disorder:
                    return {
                        "uniform_overlap": record["uniform_overlap"],
                        "rhohat_path_average": candidate[
                            "rhohat_path_average"
                        ],
                        "endpoint_overlap_at_minus_lambda": candidate[
                            "endpoint_overlap_at_minus_lambda"
                        ],
                    }
    return None


def run(args: argparse.Namespace) -> dict:
    if len(args.orders) != len(args.sample_counts):
        raise ValueError("orders and sample-counts must have equal lengths")
    if any(order % 2 for order in args.orders):
        raise ValueError("this held-out audit uses equal balanced splits")
    mp.mp.dps = args.mp_dps
    rng = np.random.default_rng(args.seed)
    spaces = {
        order // 2: exact.build_signing_space(
            order // 2, args.signing_batch_size
        )
        for order in args.orders
    }
    records = []
    for total_order, sample_count in zip(args.orders, args.sample_counts):
        child_order = total_order // 2
        for beta in args.betas:
            classes, certificate = exact.thermal_minimizer_classes(
                spaces[child_order],
                format(beta, ".12g"),
                total_order,
            )
            for left_class, left_selected in enumerate(classes):
                for right_class, right_selected in enumerate(classes):
                    left_matrix = np.asarray(
                        left_selected["representative_matrix"], dtype=np.int8
                    )
                    right_matrix = np.asarray(
                        right_selected["representative_matrix"], dtype=np.int8
                    )
                    pressure, cavity = sample_pressure_and_cavity(
                        left_matrix,
                        right_matrix,
                        beta,
                        total_order,
                        args.orientation,
                        sample_count,
                        args.batch_size,
                        rng,
                    )
                    paths = []
                    for inverse_disorder in args.lambdas:
                        path = negative_path(
                            pressure,
                            cavity,
                            inverse_disorder,
                            args.quadrature_order,
                        )
                        if left_class == 0 and right_class == 0:
                            reference = exact_reference(
                                args.exact_reference,
                                total_order,
                                beta,
                                args.orientation,
                                inverse_disorder,
                            )
                        else:
                            reference = None
                        if reference is not None:
                            path["complete_cube_reference"] = reference
                            path["sample_minus_reference"] = {
                                key: path[key] - value
                                for key, value in reference.items()
                            }
                        paths.append(path)
                    records.append(
                        {
                            "N": total_order,
                            "split": [child_order, child_order],
                            "beta": beta,
                            "raw_t": beta / math.sqrt(total_order),
                            "orientation": args.orientation,
                            "sample_count": sample_count,
                            "left_child_class": left_class,
                            "right_child_class": right_class,
                            "left_child_representative_sha256": left_selected[
                                "representative_sha256"
                            ],
                            "right_child_representative_sha256": right_selected[
                                "representative_sha256"
                            ],
                            "child_minimizer_certificate": certificate,
                            "sample_pressure_range": [
                                float(np.min(pressure)),
                                float(np.max(pressure)),
                            ],
                            "sample_pointwise_overlap_range": [
                                float(np.min(cavity)),
                                float(np.max(cavity)),
                            ],
                            "paths": paths,
                        }
                    )
                    print(
                        f"N={total_order} beta={beta:g} "
                        f"classes={left_class},{right_class} "
                        f"samples={sample_count} "
                        + " ".join(
                            f"rhohat({path['lambda']:g})="
                            f"{path['rhohat_path_average']:.8g}+-"
                            f"{path['rhohat_path_standard_error']:.2g}"
                            for path in paths
                        ),
                        flush=True,
                    )
    return {
        "schema": "actual-child-raw-negative-overlap-sample-v1",
        "classification": (
            "exhaustive thermal-child selection and exact finite Gibbs sums "
            "on reproducibly sampled fair bridges; numerical evidence only"
        ),
        "scope": {
            "orders": args.orders,
            "sample_counts": args.sample_counts,
            "betas": args.betas,
            "lambdas": args.lambdas,
            "orientation": args.orientation,
            "seed": args.seed,
            "quadrature_order": args.quadrature_order,
            "mp_dps_for_child_selection": args.mp_dps,
        },
        "exclusions": [
            "no conference or Paley surrogate",
            "no ground-state child surrogate",
            "no asymptotic conclusion from Monte Carlo values",
        ],
        "records": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--orders", type=int, nargs="+", default=[8, 10, 12, 14]
    )
    parser.add_argument(
        "--sample-counts",
        type=int,
        nargs="+",
        default=[100000, 100000, 50000, 10000],
    )
    parser.add_argument(
        "--betas", type=float, nargs="+", default=[1.0, 2.0, 4.0]
    )
    parser.add_argument("--lambdas", type=float, nargs="+", default=[1.0])
    parser.add_argument("--orientation", type=int, default=-1)
    parser.add_argument("--seed", type=int, default=20260818)
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--quadrature-order", type=int, default=32)
    parser.add_argument("--mp-dps", type=int, default=80)
    parser.add_argument("--signing-batch-size", type=int, default=8192)
    parser.add_argument(
        "--exact-reference",
        type=Path,
        default=ROOT
        / "computations/results/actual_child_negative_overlap_exact.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT
        / "computations/results/actual_child_negative_overlap_sample.json",
    )
    args = parser.parse_args()
    payload = run(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
