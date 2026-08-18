#!/usr/bin/env python3
"""Exact bridge-cube audit of the actual parent response norm in SH.0c.

Children are selected by complete signing enumeration.  Every bridge is
enumerated.  Gibbs response norms are floating-point evaluations of exact
finite sums; this is exploratory finite evidence, not an asymptotic result.
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


def weighted_quantiles(values: np.ndarray, quantiles: tuple[float, ...]) -> dict:
    return {
        format(q, ".6g"): float(np.quantile(values, q, method="higher"))
        for q in quantiles
    }


def response_norms(
    left: np.ndarray,
    right: np.ndarray,
    beta: float,
    total_n: int,
    epsilon: int,
    batch_size: int,
) -> np.ndarray:
    """Return ||E[tau X Y^T]||_F for every bridge mask."""
    m, n = len(left), len(right)
    d = m * n
    t = beta / math.sqrt(total_n)
    x = exact.projective_spins(m).astype(np.int16)
    y = exact.projective_spins(n).astype(np.int16)
    ex = exact.energies_for_matrix(left, x)
    ey = exact.energies_for_matrix(right, y)
    patterns = []
    weights = []
    for xi, exi in zip(x, ex):
        for yj, eyj in zip(y, ey):
            patterns.append((xi[:, None] * yj[None, :]).reshape(-1))
            weights.append(math.cosh(t * float(exi + epsilon * eyj)))
    q = np.asarray(patterns, dtype=np.float64)
    weight = np.asarray(weights, dtype=np.float64)
    total = 1 << d
    result = np.empty(total, dtype=np.float64)
    positions = np.arange(d, dtype=np.uint64)
    for start in range(0, total, batch_size):
        masks = np.arange(start, min(total, start + batch_size), dtype=np.uint64)
        b = 1.0 - 2.0 * (((masks[:, None] >> positions) & 1).astype(np.float64))
        field = b @ q.T
        denominator = np.cosh(t * field) @ weight
        numerator = (np.sinh(t * field) * weight[None, :]) @ q
        mean_q = numerator / denominator[:, None]
        result[start : start + len(masks)] = np.linalg.norm(mean_q, axis=1)
    return result


def run(args: argparse.Namespace) -> dict:
    mp.mp.dps = args.mp_dps
    spaces = {
        n: exact.build_signing_space(n)
        for n in range(2, args.max_total_order // 2 + 2)
    }
    cache = {}
    records = []
    for total_n in range(args.min_total_order, args.max_total_order + 1):
        m = total_n // 2
        n = total_n - m
        for beta in args.betas:
            classes = []
            for order in (m, n):
                key = (order, beta, total_n)
                if key not in cache:
                    cache[key] = exact.thermal_minimizer_classes(
                        spaces[order], format(beta, ".12g"), total_n
                    )[0]
                classes.append(cache[key])
            for left_class, right_class in itertools.product(*classes):
                left = np.asarray(left_class["representative_matrix"], dtype=np.int8)
                right = np.asarray(right_class["representative_matrix"], dtype=np.int8)
                for epsilon in (-1, 1):
                    norms = response_norms(
                        left,
                        right,
                        beta,
                        total_n,
                        epsilon,
                        args.batch_size,
                    )
                    normalized = norms / total_n
                    records.append(
                        {
                            "N": total_n,
                            "split": [m, n],
                            "beta": beta,
                            "epsilon": epsilon,
                            "left_class": left_class["class_id"],
                            "right_class": right_class["class_id"],
                            "bridge_cube_size": len(norms),
                            "response_norm_mean": float(np.mean(norms)),
                            "response_norm_maximum": float(np.max(norms)),
                            "response_norm_over_N_mean": float(np.mean(normalized)),
                            "response_norm_over_N_maximum": float(np.max(normalized)),
                            "response_norm_over_N_quantiles": weighted_quantiles(
                                normalized, (0.5, 0.9, 0.99, 0.999)
                            ),
                            "mass_response_over_N_above": {
                                format(threshold, ".6g"): float(
                                    np.mean(normalized > threshold)
                                )
                                for threshold in (0.05, 0.1, 0.2, 0.4)
                            },
                            "classification": (
                                "complete bridge enumeration; numerical Gibbs "
                                "response from exact finite spin sums"
                            ),
                        }
                    )
                    print(
                        f"N={total_n} beta={beta:g} eps={epsilon:+d} "
                        f"mean/N={np.mean(normalized):.6g} "
                        f"q99/N={np.quantile(normalized, .99):.6g} "
                        f"max/N={np.max(normalized):.6g}",
                        flush=True,
                    )
    return {
        "schema": "actual-child-parent-response-norm-exact-v1",
        "classification": (
            "complete child-signing and bridge enumeration at the declared "
            "orders; floating-point evaluation of finite Gibbs sums"
        ),
        "scope": {
            "orders": [args.min_total_order, args.max_total_order],
            "splits": "balanced only",
            "betas": args.betas,
            "orientations": [-1, 1],
            "mp_dps_for_child_selection": args.mp_dps,
        },
        "records": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-total-order", type=int, default=4)
    parser.add_argument("--max-total-order", type=int, default=8)
    parser.add_argument("--betas", type=float, nargs="+", default=[0.5, 1, 2, 4])
    parser.add_argument("--batch-size", type=int, default=65536)
    parser.add_argument("--mp-dps", type=int, default=60)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "computations/results/actual_child_response_norm_exact.json",
    )
    args = parser.parse_args()
    result = run(args)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()

