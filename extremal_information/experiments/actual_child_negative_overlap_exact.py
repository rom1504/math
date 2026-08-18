#!/usr/bin/env python3
"""Exact finite bridge-cube audit of the raw negative-tilt cavity overlap.

The children are selected by exhaustive signing/histogram enumeration at the
contracted temperature ``beta/sqrt(N)``.  For every selected child pair and
relative orientation, every bridge is enumerated.  If ``L`` is the resulting
bridge pressure and ``B^e`` flips bridge edge ``e``, the exact cavity identity

    r_e(B_-e)^2 = tanh((L(B)-L(B^e))/2)^2 / tanh(beta/sqrt(N))^2

computes the requested observable without reconstructing a surrogate Gibbs
law.  Gauss--Legendre quadrature is used only for the one-dimensional tilt
integral.  Thus child and bridge enumeration are complete, while the
transcendental evaluations and quadrature are numerical rather than interval
certificates.
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


def cavity_square_density(pressure: np.ndarray, raw_t: float) -> np.ndarray:
    """Return ``(mn)^-1 sum_e r_e(B_-e)^2`` for every bridge ``B``."""

    length = len(pressure)
    edge_count = int(round(math.log2(length)))
    if 1 << edge_count != length:
        raise ValueError("pressure array is not indexed by a Boolean cube")
    denominator = math.tanh(raw_t) ** 2
    if denominator == 0:
        raise ValueError("beta must be nonzero")
    indices = np.arange(length, dtype=np.uint64)
    density = np.zeros(length, dtype=np.float64)
    for edge in range(edge_count):
        flipped = indices ^ np.uint64(1 << edge)
        half_difference = 0.5 * (pressure - pressure[flipped])
        density += np.tanh(half_difference) ** 2 / denominator
    density /= edge_count

    # A cavity response is in [-1,1] and is unchanged by flipping its deleted
    # edge.  These checks catch both a missing half and an indexing mismatch.
    if float(np.max(density)) > 1.0 + 2e-9:
        raise AssertionError(float(np.max(density)))
    for edge in range(edge_count):
        half_difference = 0.5 * (
            pressure - pressure[indices ^ np.uint64(1 << edge)]
        )
        edge_square = np.tanh(half_difference) ** 2 / denominator
        if not np.allclose(
            edge_square,
            edge_square[indices ^ np.uint64(1 << edge)],
            rtol=0.0,
            atol=2e-12,
        ):
            raise AssertionError(f"edge-pair cavity mismatch at {edge}")
    return density


def direct_cavity_audit(
    left: np.ndarray,
    right: np.ndarray,
    beta: float,
    total_order: int,
    orientation: int,
    pressure_density: np.ndarray,
) -> dict:
    """Check selected flip-derived values against direct finite Gibbs sums."""

    rows, columns = len(left), len(right)
    edge_count = rows * columns
    raw_t = beta / math.sqrt(total_order)
    channel = math.tanh(raw_t)
    x = exact.projective_spins(rows).astype(np.int16)
    y = exact.projective_spins(columns).astype(np.int16)
    ex = exact.energies_for_matrix(left, x)
    ey = exact.energies_for_matrix(right, y)
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
    cube_size = 1 << edge_count
    masks = sorted(set([0, cube_size // 3, cube_size - 1]))
    errors = []
    direct_values = []
    positions = np.arange(edge_count, dtype=np.uint64)
    for mask in masks:
        bridge = 1.0 - 2.0 * (
            ((np.uint64(mask) >> positions) & 1).astype(np.float64)
        )
        field = rank_one @ bridge
        denominator = float(
            np.dot(internal, np.cosh(raw_t * field))
        )
        full_response = (
            (internal * np.sinh(raw_t * field)) @ rank_one
        ) / denominator
        cavity_response = (
            full_response - channel * bridge
        ) / (1.0 - channel * bridge * full_response)
        direct = float(np.mean(cavity_response * cavity_response))
        direct_values.append(direct)
        errors.append(abs(direct - float(pressure_density[mask])))
    maximum_error = max(errors)
    if maximum_error > 2e-9:
        raise AssertionError((masks, direct_values, maximum_error))
    return {
        "checked_bridge_masks": masks,
        "direct_normalized_cavity_square_values": direct_values,
        "maximum_flip_vs_direct_absolute_error": maximum_error,
    }


def tilted_mean(
    values: np.ndarray, pressure: np.ndarray, tilt: float
) -> tuple[float, float]:
    exponent = tilt * pressure
    exponent -= float(np.max(exponent))
    weight = np.exp(exponent)
    total = float(np.sum(weight))
    mean = float(np.dot(weight, values) / total)
    effective_sample_fraction = total * total / (
        len(weight) * float(np.dot(weight, weight))
    )
    return mean, effective_sample_fraction


def path_record(
    pressure: np.ndarray,
    cavity_density: np.ndarray,
    inverse_disorder: float,
    quadrature_order: int,
) -> dict:
    nodes, weights = np.polynomial.legendre.leggauss(quadrature_order)
    integral = 0.0
    minimum_effective_fraction = 1.0
    for node, weight in zip(nodes, weights):
        tilt = 0.5 * inverse_disorder * (float(node) - 1.0)
        mean, effective_fraction = tilted_mean(
            cavity_density, pressure, tilt
        )
        integral += 0.5 * float(weight) * mean
        minimum_effective_fraction = min(
            minimum_effective_fraction, effective_fraction
        )
    endpoint, endpoint_effective_fraction = tilted_mean(
        cavity_density, pressure, -inverse_disorder
    )
    midpoint, midpoint_effective_fraction = tilted_mean(
        cavity_density, pressure, -0.5 * inverse_disorder
    )
    return {
        "lambda": inverse_disorder,
        "rhohat_path_average": integral,
        "endpoint_overlap_at_minus_lambda": endpoint,
        "midpoint_overlap_at_minus_lambda_over_2": midpoint,
        "endpoint_importance_effective_fraction": endpoint_effective_fraction,
        "midpoint_importance_effective_fraction": midpoint_effective_fraction,
        "minimum_quadrature_importance_effective_fraction": (
            minimum_effective_fraction
        ),
        "quadrature_order": quadrature_order,
    }


def run(args: argparse.Namespace) -> dict:
    mp.mp.dps = args.mp_dps
    order_plans = []
    for total_order in range(args.min_total_order, args.max_total_order + 1):
        left_order = (
            args.left_order
            if args.left_order is not None
            else total_order // 2
        )
        right_order = total_order - left_order
        if not (2 <= left_order <= right_order):
            raise ValueError(
                f"invalid split {left_order}+{right_order} for N={total_order}"
            )
        order_plans.append((total_order, left_order, right_order))
    largest_child = max(max(left, right) for _, left, right in order_plans)
    spaces = {
        order: exact.build_signing_space(order, args.signing_batch_size)
        for order in range(2, largest_child + 1)
    }
    minimizer_cache = {}
    child_certificates = []
    records = []
    for total_order, left_order, right_order in order_plans:
        for beta in args.betas:
            beta_text = format(beta, ".12g")
            class_lists = []
            for order in (left_order, right_order):
                key = (order, beta_text, total_order)
                if key not in minimizer_cache:
                    minimizer_cache[key] = exact.thermal_minimizer_classes(
                        spaces[order], beta_text, total_order
                    )
                    child_certificates.append(minimizer_cache[key][1])
                class_lists.append(minimizer_cache[key][0])
            for left_class, right_class in itertools.product(*class_lists):
                left = np.asarray(
                    left_class["representative_matrix"], dtype=np.int8
                )
                right = np.asarray(
                    right_class["representative_matrix"], dtype=np.int8
                )
                for orientation in args.orientations:
                    pressure, pressure_audit = exact.bridge_pressures(
                        left, right, beta, total_order, orientation
                    )
                    raw_t = beta / math.sqrt(total_order)
                    density = cavity_square_density(pressure, raw_t)
                    direct_audit = direct_cavity_audit(
                        left,
                        right,
                        beta,
                        total_order,
                        orientation,
                        density,
                    )
                    uniform_mean, _ = tilted_mean(density, pressure, 0.0)
                    records.append(
                        {
                            "N": total_order,
                            "split": [left_order, right_order],
                            "beta": beta,
                            "raw_t": raw_t,
                            "orientation": orientation,
                            "left_child_class": left_class["class_id"],
                            "right_child_class": right_class["class_id"],
                            "left_child_sha256": left_class[
                                "representative_sha256"
                            ],
                            "right_child_sha256": right_class[
                                "representative_sha256"
                            ],
                            "uniform_overlap": uniform_mean,
                            "pointwise_overlap_minimum": float(
                                np.min(density)
                            ),
                            "pointwise_overlap_maximum": float(
                                np.max(density)
                            ),
                            "paths": [
                                path_record(
                                    pressure,
                                    density,
                                    inverse_disorder,
                                    args.quadrature_order,
                                )
                                for inverse_disorder in args.lambdas
                            ],
                            "bridge_pressure_audit": pressure_audit,
                            "cavity_direct_audit": direct_audit,
                        }
                    )
                    print(
                        f"N={total_order} beta={beta:g} eps={orientation:+d} "
                        f"rho0={uniform_mean:.8g} "
                        + " ".join(
                            f"rhohat({path['lambda']:g})="
                            f"{path['rhohat_path_average']:.8g}"
                            for path in records[-1]["paths"]
                        ),
                        flush=True,
                    )
    return {
        "schema": "actual-child-raw-negative-overlap-exact-v1",
        "classification": (
            "complete child-signing and bridge enumeration; numerical "
            "transcendental evaluation and Gauss-Legendre quadrature"
        ),
        "normalization": (
            "rhohat=(lambda*mn)^-1 int_-lambda^0 "
            "E_Pi_s sum_e r_e(B_-e)^2 ds"
        ),
        "scope": {
            "orders": [args.min_total_order, args.max_total_order],
            "splits": (
                "balanced floor(N/2)+ceil(N/2)"
                if args.left_order is None
                else f"fixed left order {args.left_order}"
            ),
            "betas": args.betas,
            "lambdas": args.lambdas,
            "orientations": args.orientations,
            "quadrature_order": args.quadrature_order,
            "mp_dps_for_child_selection": args.mp_dps,
        },
        "exclusions": [
            "no conference or Paley surrogate",
            "no ground-state child surrogate",
            "no bridge sampling",
            "no asymptotic inference from finite values",
        ],
        "child_minimizer_certificates": child_certificates,
        "records": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-total-order", type=int, default=4)
    parser.add_argument("--max-total-order", type=int, default=9)
    parser.add_argument(
        "--left-order",
        type=int,
        default=None,
        help="override the balanced left order (useful for exact 3+7 at N=10)",
    )
    parser.add_argument(
        "--betas", type=float, nargs="+", default=[1.0, 2.0, 4.0]
    )
    parser.add_argument(
        "--lambdas",
        type=float,
        nargs="+",
        default=[0.5, 1.0, 2.0, 4.0, 5.382104195764755],
    )
    parser.add_argument(
        "--orientations", type=int, nargs="+", default=[-1, 1]
    )
    parser.add_argument("--quadrature-order", type=int, default=32)
    parser.add_argument("--mp-dps", type=int, default=80)
    parser.add_argument("--signing-batch-size", type=int, default=8192)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT
        / "computations/results/actual_child_negative_overlap_exact.json",
    )
    args = parser.parse_args()
    result = run(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
