#!/usr/bin/env python3
"""Finite verifier for the actual-child row-ANOVA tangent formulas.

This is a normalization and algebra check for Theorems RA.1--RA.2 in
``drafts/actual_child_row_anova_infinitesimal.md``.  It enumerates the child
signings and bridge cube completely at one small order.  Transcendental
quantities and coordinate-product minimization are numerical; no asymptotic
claim is inferred from the output.
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
import actual_child_row_product_shadow as shadow  # noqa: E402


def row_anova(pressure: np.ndarray, rows: int, columns: int) -> dict:
    indices = np.arange(len(pressure), dtype=np.uint64)
    mean = float(np.mean(pressure))
    centered = pressure - mean
    additive = np.zeros_like(pressure)
    row_norms = []
    mask = (1 << columns) - 1
    for row in range(rows):
        code = ((indices >> (row * columns)) & mask).astype(np.int64)
        sums = np.bincount(code, weights=pressure, minlength=1 << columns)
        counts = np.bincount(code, minlength=1 << columns)
        component = sums / counts - mean
        additive += component[code]
        row_norms.append(float(np.mean(component * component)))
    total = float(np.mean(centered * centered))
    additive_variance = float(np.mean(additive * additive))
    cross = float(np.mean((centered - additive) ** 2))
    return {
        "total_variance": total,
        "additive_variance": additive_variance,
        "cross_variance": cross,
        "row_component_variances": row_norms,
        "orthogonality_residual": total - additive_variance - cross,
    }


def mixed_row_response(pressure: np.ndarray, rows: int, columns: int) -> dict:
    """Compute J_2 by exact double-centering of every pair of row axes."""
    tensor = shadow.pressure_tensor(pressure, rows, 1 << columns)
    pair_masses = []
    for i in range(rows):
        for k in range(i + 1, rows):
            moved = np.moveaxis(tensor, (i, k), (0, 1))
            matrix_stack = moved.reshape(moved.shape[0], moved.shape[1], -1)
            double_centered = (
                matrix_stack
                - np.mean(matrix_stack, axis=0, keepdims=True)
                - np.mean(matrix_stack, axis=1, keepdims=True)
                + np.mean(matrix_stack, axis=(0, 1), keepdims=True)
            )
            pair_masses.append(
                {
                    "rows": [i, k],
                    "quarter_mean_squared_mixed_difference": float(
                        np.mean(double_centered * double_centered)
                    ),
                }
            )
    return {
        "pair_masses": pair_masses,
        "J2": sum(
            row["quarter_mean_squared_mixed_difference"] for row in pair_masses
        ),
    }


def generalized_bridge_pressures(
    left: np.ndarray,
    right: np.ndarray,
    internal_t: float,
    bridge_u: float,
    epsilon: int,
) -> np.ndarray:
    """Exact bridge cube with separate internal and bridge amplitudes."""
    m, n = len(left), len(right)
    d = m * n
    length = 1 << d
    x = exact.projective_spins(m).astype(np.int16)
    y = exact.projective_spins(n).astype(np.int16)
    ex = exact.energies_for_matrix(left, x)
    ey = exact.energies_for_matrix(right, y)
    weights = np.zeros(length, dtype=np.longdouble)
    for xi, exi in zip(x, ex):
        for yj, eyj in zip(y, ey):
            signs = (xi[:, None] * yj[None, :]).reshape(-1)
            pattern = 0
            for bit, sign in enumerate(signs):
                if sign < 0:
                    pattern |= 1 << bit
            weights[pattern] += np.cosh(
                np.longdouble(internal_t)
                * np.longdouble(exi + epsilon * eyj)
            )
    indices = np.arange(length, dtype=np.uint64)
    popcount = np.bitwise_count(indices).astype(np.int16)
    radial = np.cosh(
        np.longdouble(bridge_u) * (d - 2 * popcount).astype(np.longdouble)
    )
    zbar = exact.xor_convolution(weights, radial) / np.longdouble(len(x) * len(y))
    if np.min(zbar) <= 0:
        raise FloatingPointError("nonpositive generalized partition value")
    return np.log(zbar).astype(np.float64)


def overlap_cross_coefficient(
    left: np.ndarray,
    right: np.ndarray,
    internal_t: float,
    epsilon: int,
) -> dict:
    x = exact.projective_spins(len(left)).astype(np.float64)
    y = exact.projective_spins(len(right)).astype(np.float64)
    ex = exact.energies_for_matrix(left, x.astype(np.int16))
    ey = exact.energies_for_matrix(right, y.astype(np.int16))
    weight = np.cosh(
        internal_t * (ex[:, None] + epsilon * ey[None, :])
    )
    weight /= np.sum(weight)
    cross = 0.0
    same_column_lower = 0.0
    for i in range(len(left)):
        for k in range(i + 1, len(left)):
            cx = x[:, i] * x[:, k]
            for j in range(len(right)):
                for ell in range(len(right)):
                    cy = y[:, j] * y[:, ell]
                    gamma = float(np.einsum("a,b,ab->", cx, cy, weight))
                    cross += gamma * gamma
                    if j == ell:
                        same_column_lower += gamma * gamma
    additive = 0.0
    for i in range(len(left)):
        for j in range(len(right)):
            for ell in range(j + 1, len(right)):
                cy = y[:, j] * y[:, ell]
                gamma = float(np.einsum("ab,b->", weight, cy))
                additive += gamma * gamma
    return {
        "K_cross": cross,
        "same_column_lower_bound": same_column_lower,
        "K_additive": additive,
    }


def escort_metrics(pressure: np.ndarray, rows: int, columns: int, lam: float) -> dict:
    shifted = -lam * pressure
    shifted -= np.max(shifted)
    q = np.exp(shifted)
    q /= np.sum(q)
    uniform_mean = float(np.mean(pressure))
    soft = exact.negative_moment_soft_pressure(pressure, lam)
    tc = exact.block_information_metrics(q, rows, columns)["row_total_correlation"]
    product = shadow.optimize_product_shadow(
        pressure,
        rows,
        columns,
        lam,
        random_starts=0,
        seed=0,
        tolerance=1e-15,
        max_sweeps=10000,
    )
    return {
        "lambda": lam,
        "G": uniform_mean - soft,
        "row_product_gain": product["rigorous_candidate_row_product_gain"],
        "candidate_I_left": product["candidate_reverse_projection_upper_bound"],
        "row_total_correlation": tc,
    }


def run(args: argparse.Namespace) -> dict:
    mp.mp.dps = 80
    m = args.left_order
    n = args.total_order - m
    spaces = {k: exact.build_signing_space(k) for k in {m, n}}
    classes = {}
    for k in (m, n):
        classes[k] = exact.thermal_minimizer_classes(
            spaces[k], format(args.beta, ".12g"), args.total_order
        )[0]
    left = np.asarray(classes[m][0]["representative_matrix"], dtype=np.int8)
    right = np.asarray(classes[n][0]["representative_matrix"], dtype=np.int8)
    t = args.beta / math.sqrt(args.total_order)
    pressure = generalized_bridge_pressures(left, right, t, t, args.epsilon)
    anova = row_anova(pressure, m, n)
    mixed = mixed_row_response(pressure, m, n)
    overlap = overlap_cross_coefficient(left, right, t, args.epsilon)
    amplitudes = []
    for u in args.amplitudes:
        value = row_anova(
            generalized_bridge_pressures(left, right, t, u, args.epsilon), m, n
        )
        amplitudes.append(
            {
                "u": u,
                "cross_variance_over_u4": value["cross_variance"] / u**4,
                "additive_variance_over_u4": value["additive_variance"] / u**4,
            }
        )
    laws = [escort_metrics(pressure, m, n, lam) for lam in args.lambdas]
    for law in laws:
        lam = law["lambda"]
        law["G_over_lambda_half_variance"] = law["G"] / (
            0.5 * lam * anova["total_variance"]
        )
        law["product_gain_over_lambda_half_additive"] = law[
            "row_product_gain"
        ] / (0.5 * lam * anova["additive_variance"])
        law["I_left_over_lambda2_half_cross"] = law["candidate_I_left"] / (
            0.5 * lam * lam * anova["cross_variance"]
        )
        law["TC_over_lambda2_half_cross"] = law["row_total_correlation"] / (
            0.5 * lam * lam * anova["cross_variance"]
        )
    return {
        "schema": "actual-child-row-anova-infinitesimal-verifier-v1",
        "classification": (
            "complete finite enumeration; numerical transcendental evaluation; "
            "the tiny-lambda row-product optimizer is obtained by coordinate "
            "best responses from the uniform law"
        ),
        "N": args.total_order,
        "split": [m, n],
        "beta": args.beta,
        "raw_t": t,
        "epsilon": args.epsilon,
        "left_sha": exact.matrix_sha(left),
        "right_sha": exact.matrix_sha(right),
        "physical_row_anova": anova,
        "physical_mixed_row_response": mixed,
        "mixed_response_lower_bound_residual": (
            mixed["J2"] - anova["cross_variance"]
        ),
        "mixed_response_upper_bound_residual": (
            math.comb(m, 2) * anova["cross_variance"] - mixed["J2"]
        ),
        "overlap_coefficients": overlap,
        "bridge_amplitude_checks": amplitudes,
        "disorder_temperature_checks": laws,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--total-order", type=int, default=6)
    parser.add_argument("--left-order", type=int, default=3)
    parser.add_argument("--beta", type=float, default=1.0)
    parser.add_argument("--epsilon", type=int, choices=(-1, 1), default=1)
    parser.add_argument(
        "--amplitudes", type=float, nargs="+", default=[0.1, 0.05, 0.025, 0.0125]
    )
    parser.add_argument(
        "--lambdas", type=float, nargs="+", default=[0.1, 0.05, 0.025, 0.0125]
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "computations/results/actual_child_row_anova_verify.json",
    )
    args = parser.parse_args()
    result = run(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
