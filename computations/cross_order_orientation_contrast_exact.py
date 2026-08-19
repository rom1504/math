#!/usr/bin/env python3
"""Exact small-order audit of the free relative-orientation gain.

For each split and beta, this program selects representatives of every exact
thermal-minimizer class at the *child's own* scale beta/sqrt(child_order),
enumerates every sign bridge in both relative orientations, and evaluates the
two direct cross-order certificates

    (E_B L_+ + E_B L_-)/2 - P_m - P_n,
    E_B min(L_+, L_-) - P_m - P_n.

Signing and bridge cubes and all integer energies are enumerated exactly.
Only transcendental pressure evaluations are numerical; child histogram
comparisons use mpmath at the requested precision.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import mpmath as mp
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from extremal_information.experiments.actual_child_bridge_law_exact import (
    bridge_pressures,
    build_signing_space,
    thermal_minimizer_classes,
)


def joint_soft_min(
    minus: np.ndarray,
    plus: np.ndarray,
    lam: float,
) -> float:
    pivot = min(float(np.min(minus)), float(np.min(plus)))
    moment = 0.5 * (
        float(np.mean(np.exp(-lam * (minus - pivot))))
        + float(np.mean(np.exp(-lam * (plus - pivot))))
    )
    return pivot - math.log(moment) / lam


def audit(
    min_total_n: int,
    max_total_n: int,
    betas: list[float],
    mp_dps: int,
    balanced_only: bool,
) -> dict:
    mp.mp.dps = mp_dps
    needed_orders = set()
    for total_n in range(min_total_n, max_total_n + 1):
        splits = (
            [(total_n // 2, total_n - total_n // 2)]
            if balanced_only
            else [(m, total_n - m) for m in range(2, total_n // 2 + 1)]
        )
        for m, n in splits:
            if m >= 2:
                needed_orders.update((m, n))
    spaces = {k: build_signing_space(k) for k in sorted(needed_orders)}
    cache: dict[tuple[int, str], tuple[list[dict], dict]] = {}

    def minimizers(k: int, beta_text: str) -> tuple[list[dict], dict]:
        key = (k, beta_text)
        if key not in cache:
            cache[key] = thermal_minimizer_classes(
                spaces[k], beta_text, k
            )
        return cache[key]

    records: list[dict] = []
    for total_n in range(min_total_n, max_total_n + 1):
        splits = (
            [(total_n // 2, total_n - total_n // 2)]
            if balanced_only
            else [(m, total_n - m) for m in range(2, total_n // 2 + 1)]
        )
        for m, n in splits:
            n = total_n - m
            for beta in betas:
                beta_text = format(beta, ".12g")
                left_classes, left_certificate = minimizers(m, beta_text)
                right_classes, right_certificate = minimizers(n, beta_text)
                child_target = float(left_certificate["mp_optimum"]) + float(
                    right_certificate["mp_optimum"]
                )
                for left_class in left_classes:
                    left = np.asarray(
                        left_class["representative_matrix"], dtype=np.int8
                    )
                    for right_class in right_classes:
                        right = np.asarray(
                            right_class["representative_matrix"], dtype=np.int8
                        )
                        minus, minus_audit = bridge_pressures(
                            left, right, beta, total_n, -1
                        )
                        plus, plus_audit = bridge_pressures(
                            left, right, beta, total_n, +1
                        )
                        contrast = 0.5 * float(np.mean(np.abs(plus - minus)))
                        uniform_defect = (
                            0.5 * float(np.mean(plus) + np.mean(minus))
                            - child_target
                        )
                        orientation_selected_defect = (
                            float(np.mean(np.minimum(plus, minus)))
                            - child_target
                        )
                        best_bridge_defect = (
                            min(float(np.min(plus)), float(np.min(minus)))
                            - child_target
                        )
                        basin_thresholds = (0.0, 0.05, 0.1, 0.25, 0.5, 1.0)
                        joint_basin_mass = {
                            format(threshold, ".12g"): 0.5
                            * (
                                float(np.mean(minus <= child_target + threshold))
                                + float(np.mean(plus <= child_target + threshold))
                            )
                            for threshold in basin_thresholds
                        }
                        fractional_certificates = {
                            format(lam, ".12g"): (
                                joint_soft_min(minus, plus, lam) - child_target
                            )
                            for lam in (0.25, 0.5, 1.0, 2.0, 4.0, 8.0)
                        }
                        identity_error = abs(
                            orientation_selected_defect
                            - (uniform_defect - contrast)
                        )
                        if identity_error > 2e-10:
                            raise AssertionError(identity_error)
                        records.append(
                            {
                                "N": total_n,
                                "split": [m, n],
                                "beta": beta,
                                "left_child_class": left_class["class_id"],
                                "right_child_class": right_class["class_id"],
                                "left_child_sha256": left_class[
                                    "representative_sha256"
                                ],
                                "right_child_sha256": right_class[
                                    "representative_sha256"
                                ],
                                "child_own_scale_pressure_sum": child_target,
                                "uniform_joint_defect": uniform_defect,
                                "mean_half_absolute_orientation_contrast": contrast,
                                "orientation_selected_mean_defect": (
                                    orientation_selected_defect
                                ),
                                "best_bridge_and_orientation_defect": (
                                    best_bridge_defect
                                ),
                                "joint_near_target_basin_mass": joint_basin_mass,
                                "fractional_soft_min_defect": fractional_certificates,
                                "contrast_fraction_of_uniform_defect": (
                                    contrast / uniform_defect
                                    if uniform_defect > 0
                                    else None
                                ),
                                "selected_identity_absolute_error": identity_error,
                                "minus_bridge_audit": minus_audit,
                                "plus_bridge_audit": plus_audit,
                            }
                        )
                        print(
                            f"N={total_n} {m}+{n} beta={beta:g} "
                            f"classes={left_class['class_id']},"
                            f"{right_class['class_id']} "
                            f"uniform={uniform_defect:.6g} "
                            f"contrast={contrast:.6g} "
                            f"best={best_bridge_defect:.6g}",
                            flush=True,
                        )
    return {
        "schema": "cross-order-orientation-contrast-exact-v1",
        "classification": (
            "exact signing/bridge enumeration and integer energies; "
            "high-precision numerical child selection and floating-point "
            "transcendental pressure evaluation"
        ),
        "normalization": (
            "P_k(beta)=min_A log(2^-k sum_x cosh(beta H_A(x)/sqrt(k))); "
            "parent bridge scale beta/sqrt(m+n)"
        ),
        "scope": {
            "min_total_n": min_total_n,
            "max_total_n": max_total_n,
            "betas": betas,
            "mp_dps": mp_dps,
            "balanced_only": balanced_only,
            "child_selection": "exact thermal minimizers at child-own scale",
        },
        "records": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--min-total-n", type=int, default=4)
    parser.add_argument("--max-total-n", type=int, default=9)
    parser.add_argument("--balanced-only", action="store_true")
    parser.add_argument(
        "--betas", type=float, nargs="+", default=[0.25, 0.5, 1.0, 2.0, 4.0]
    )
    parser.add_argument("--mp-dps", type=int, default=80)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "computations/results/cross_order_orientation_contrast_exact.json"
        ),
    )
    args = parser.parse_args()
    result = audit(
        args.min_total_n,
        args.max_total_n,
        args.betas,
        args.mp_dps,
        args.balanced_only,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
