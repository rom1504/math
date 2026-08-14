#!/usr/bin/env python3
"""Cheap held-out full two-block basin probe at k=4,5,6.

The certified audit constructor verifies a Hadamard identity by dense matrix
multiplication.  That check is useful at small order but cubic at k=6.  This
runner constructs the identical matrices directly from their defining Walsh
formulas and reuses the audited ``full_sample`` routine.  It first checks the
fast constructor entry-for-entry at k=2 and k=3.

All random inputs are uniform after quotienting the harmless global sign by
fixing their first coordinate to +1.  The output labels Monte Carlo evidence
as such; it is not an optimality or asymptotic certificate.
"""

from __future__ import annotations

import argparse
import gc
import json
import math
import time
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, Tuple

import numpy as np

from audit_regular_hadamard_walsh_basins import (
    full_obstruction_matrix,
    full_sample,
    parity_table,
    sylvester,
)


def fast_full_obstruction_matrix(k: int) -> np.ndarray:
    """Construct the audited two-block matrix without a cubic identity check."""
    s = 1 << k
    m = s * s
    values = np.arange(m, dtype=np.int64)
    a = values & (s - 1)
    b = values >> k
    parity = parity_table(s)
    f = (1 - 2 * parity[np.bitwise_and(a, b)]).astype(np.int32)
    regular_hadamard = f[np.bitwise_xor(values[:, None], values[None, :])]
    signing = regular_hadamard - np.eye(m, dtype=np.int32)
    laplacian_edge = np.array([[1, -1], [-1, 1]], dtype=np.int32)
    bridge = np.kron(sylvester(m // 2), laplacian_edge).astype(np.int32)
    return np.block(
        [[signing, bridge], [bridge.T, -signing]]
    ).astype(np.int32)


def parse_sample_spec(specification: str) -> Iterable[Tuple[int, int]]:
    for item in specification.split(","):
        k_text, samples_text = item.split(":", maxsplit=1)
        k, samples = int(k_text), int(samples_text)
        if k < 2 or k > 6 or samples <= 0:
            raise ValueError("each item must have 2 <= k <= 6 and samples > 0")
        yield k, samples


def zero_count_upper_95(samples: int) -> float:
    """Exact one-sided 95% binomial upper endpoint after zero hits."""
    return 1.0 - 0.05 ** (1.0 / samples)


def run_order(k: int, samples: int, seed: int) -> Dict[str, object]:
    build_start = time.perf_counter()
    matrix = fast_full_obstruction_matrix(k)
    build_seconds = time.perf_counter() - build_start
    order = len(matrix)
    rng = np.random.default_rng(seed)
    thresholds = (0.0, 0.005, 0.01, 0.025, 0.05, 0.1)
    threshold_counts = Counter({threshold: 0 for threshold in thresholds})
    positive_defect_count = 0
    hard_branch_count = 0
    nonzero_kappas = []

    sample_start = time.perf_counter()
    for _ in range(samples):
        spin = rng.choice(np.array([-1, 1], dtype=np.int32), size=order)
        spin[0] = 1
        metrics = full_sample(matrix, {}, spin, None)["metrics"]
        kappa = float(metrics["kappa_over_N"])
        positive_defect_count += int(float(metrics["defect_over_N32"]) > 0.0)
        hard_branch_count += int(float(metrics["hard_branch"]) > 0.0)
        if kappa > 0.0:
            nonzero_kappas.append(kappa)
        for threshold in thresholds:
            threshold_counts[threshold] += int(kappa > threshold)
    sample_seconds = time.perf_counter() - sample_start

    threshold_report: Dict[str, object] = {}
    for threshold in thresholds:
        count = threshold_counts[threshold]
        report: Dict[str, object] = {
            "count": count,
            "empirical_probability": count / samples,
        }
        if count == 0:
            report["exact_one_sided_95_percent_upper_bound"] = (
                zero_count_upper_95(samples)
            )
        threshold_report[str(threshold)] = report

    answer: Dict[str, object] = {
        "k": k,
        "N": order,
        "samples_after_global_sign_quotient": samples,
        "seed": seed,
        "matrix_mebibytes": matrix.nbytes / (1 << 20),
        "build_seconds": build_seconds,
        "sample_seconds": sample_seconds,
        "samples_per_second": samples / sample_seconds,
        "hard_branch_count": hard_branch_count,
        "positive_defect_count": positive_defect_count,
        "positive_defect_probability": positive_defect_count / samples,
        "kappa_over_N_thresholds_strict": threshold_report,
        "nonzero_kappa_count": len(nonzero_kappas),
        "max_kappa_over_N": max(nonzero_kappas, default=0.0),
    }
    if nonzero_kappas:
        answer["nonzero_kappa_quantiles"] = {
            str(q): float(np.quantile(nonzero_kappas, q))
            for q in (0.25, 0.5, 0.75, 0.9, 0.99)
        }
    del matrix
    gc.collect()
    return answer


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--samples",
        default="4:10000,5:2000,6:250",
        help="comma-separated k:sample-count entries",
    )
    parser.add_argument("--seed", type=int, default=20260814)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "computations/results/regular_hadamard_full_basin_scaling.json"
        ),
    )
    args = parser.parse_args()

    constructor_checks: Dict[str, bool] = {}
    for k in (2, 3):
        audited, _ = full_obstruction_matrix(k)
        fast = fast_full_obstruction_matrix(k)
        constructor_checks[f"k={k}"] = bool(np.array_equal(audited, fast))
        if not constructor_checks[f"k={k}"]:
            raise AssertionError(f"fast constructor differs at k={k}")

    orders = []
    for index, (k, samples) in enumerate(parse_sample_spec(args.samples)):
        order_seed = args.seed + 1_000_003 * k + index
        print(f"running k={k}, samples={samples}, seed={order_seed}", flush=True)
        result = run_order(k, samples, order_seed)
        orders.append(result)
        print(
            "  positive defects:",
            result["positive_defect_count"],
            "max kappa/N:",
            result["max_kappa_over_N"],
            "seconds:",
            round(float(result["sample_seconds"]), 3),
            flush=True,
        )

    payload = {
        "classification": "fixed-seed Monte Carlo evidence",
        "experiment": "uniform full two-block regular-Hadamard obstruction",
        "event": "kappa_*/N exceeds each listed threshold strictly",
        "constructor_checks": constructor_checks,
        "orders": orders,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(args.output)


if __name__ == "__main__":
    main()
