#!/usr/bin/env python3
"""Spectral audit for conference children joined by a random sign bridge.

For a symmetric conference signing A of order r and an iid Rademacher
bridge B, form

    W_eps = [[A, B], [B.T, eps*A]] / sqrt(2*r),  eps in {-1,+1}.

The script measures normalized spectral moments, operator norms, the
entrywise residuals in Assumption 2.9 of arXiv:2607.10102, and two centered
mixed moments which must vanish if the deterministic and random blocks are
asymptotically scalar-free.  The two orientations for one bridge are kept as
paired observations.  All matrix power traces are exact floating-point
evaluations for the sampled finite matrices; bridge sampling is Monte Carlo.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]

SOURCES = {
    6: ("conference_double_p5.json", "conference_matrix"),
    10: ("conference_order10_gf9.json", "conference_matrix"),
    14: ("conference_double_p13.json", "conference_matrix"),
    18: ("conference_double_p17.json", "conference_matrix"),
    26: ("conference_order26_gf25.json", "conference_matrix"),
    98: ("conference_double_p97.json", "conference_matrix"),
}

DEFAULT_SAMPLES = {6: 1200, 10: 1000, 14: 800, 18: 600, 26: 400, 98: 120}


def load_conference(order: int) -> tuple[np.ndarray, str]:
    filename, key = SOURCES[order]
    path = ROOT / "computations" / "results" / filename
    payload = json.loads(path.read_text(encoding="utf-8"))
    matrix = np.asarray(payload[key], dtype=np.float64)
    target = (order - 1) * np.eye(order)
    if not np.array_equal(matrix @ matrix, target):
        raise ValueError(f"{path} is not a symmetric conference matrix")
    return matrix, str(path.relative_to(ROOT))


def summary(values: np.ndarray) -> dict:
    values = np.asarray(values, dtype=np.float64)
    count = len(values)
    return {
        "mean": float(np.mean(values)),
        "standard_error": (
            0.0 if count <= 1 else float(np.std(values, ddof=1) / math.sqrt(count))
        ),
        "q05": float(np.quantile(values, 0.05)),
        "q50": float(np.quantile(values, 0.50)),
        "q95": float(np.quantile(values, 0.95)),
        "maximum": float(np.max(values)),
    }


def theoretical_moments(order: int) -> dict:
    r = order
    return {
        "m2_exact_expectation": 1.0 - 1.0 / (2.0 * r),
        "m4_exact_expectation": (
            14 * r**3 - 14 * r**2 + 2 * r
        ) / (8 * r**3),
        "m6_exact_expectation": (
            60 * r**4 - 90 * r**3 + 34 * r**2 - 2 * r
        ) / (16 * r**4),
        "free_limit_m2": 1.0,
        "free_limit_m4": 7.0 / 4.0,
        "free_limit_m6": 15.0 / 4.0,
        "free_limit_m8": 143.0 / 16.0,
        "free_limit_outer_edge": 3.0 * math.sqrt(6.0) / 4.0,
    }


def audit_order(order: int, samples: int, seed: int) -> dict:
    a, source = load_conference(order)
    r = order
    n = 2 * r
    root_n = math.sqrt(n)
    rng = np.random.default_rng(seed)
    identity = np.eye(n)
    zero = np.zeros_like(a)

    observations = {
        epsilon: {
            "m2": [],
            "m4": [],
            "m6": [],
            "m8": [],
            "operator_norm": [],
            "mixed_DR_fourth": [],
            "mixed_centered_R2": [],
            **{
                f"power_{power}_diag_residual": []
                for power in range(1, 7)
            },
            **{
                f"power_{power}_offdiag_max": []
                for power in range(1, 7)
            },
        }
        for epsilon in (-1, 1)
    }

    for _ in range(samples):
        b = rng.choice(np.asarray([-1.0, 1.0]), size=(r, r))
        random_block = np.block([[zero, b], [b.T, zero]]) / root_n
        r2_centered = random_block @ random_block - 0.5 * identity

        for epsilon in (-1, 1):
            deterministic_block = (
                np.block([[a, zero], [zero, epsilon * a]]) / root_n
            )
            parent = deterministic_block + random_block
            powers = {1: parent}
            for power in range(2, 9):
                powers[power] = powers[power - 1] @ parent

            row = observations[epsilon]
            for power in (2, 4, 6, 8):
                row[f"m{power}"].append(float(np.trace(powers[power]) / n))
            row["operator_norm"].append(
                float(max(abs(np.linalg.eigvalsh(parent))))
            )

            dr = deterministic_block @ random_block
            row["mixed_DR_fourth"].append(
                float(np.trace(dr @ dr @ dr @ dr) / n)
            )
            row["mixed_centered_R2"].append(
                float(
                    np.trace(
                        deterministic_block
                        @ r2_centered
                        @ deterministic_block
                        @ r2_centered
                    )
                    / n
                )
            )

            for power in range(1, 7):
                matrix = powers[power]
                trace_mean = float(np.trace(matrix) / n)
                diagonal = np.diag(matrix)
                diag_residual = float(np.max(np.abs(diagonal - trace_mean)))
                off_diagonal = matrix.copy()
                np.fill_diagonal(off_diagonal, 0.0)
                offdiag_max = float(np.max(np.abs(off_diagonal)))
                row[f"power_{power}_diag_residual"].append(diag_residual)
                row[f"power_{power}_offdiag_max"].append(offdiag_max)

    scale = math.sqrt(n / math.log(n))
    output_orientations = []
    for epsilon in (-1, 1):
        row = observations[epsilon]
        result = {
            "epsilon": epsilon,
            "moments": {
                name: summary(np.asarray(row[name]))
                for name in ("m2", "m4", "m6", "m8")
            },
            "operator_norm": summary(np.asarray(row["operator_norm"])),
            "mixed_freeness_probes": {
                "tr_normalized_(D_R)^4": summary(
                    np.asarray(row["mixed_DR_fourth"])
                ),
                "tr_normalized_D_(R2-halfI)_D_(R2-halfI)": summary(
                    np.asarray(row["mixed_centered_R2"])
                ),
            },
            "assumption_2_9_power_residuals": [],
        }
        for power in range(1, 7):
            diag = np.asarray(row[f"power_{power}_diag_residual"])
            offdiag = np.asarray(row[f"power_{power}_offdiag_max"])
            result["assumption_2_9_power_residuals"].append(
                {
                    "power": power,
                    "diag_residual": summary(diag),
                    "offdiag_max": summary(offdiag),
                    "sqrt_N_over_logN_scaled_diag": summary(scale * diag),
                    "sqrt_N_over_logN_scaled_offdiag": summary(scale * offdiag),
                }
            )
        output_orientations.append(result)

    paired = {}
    for name in ("m2", "m4", "m6", "m8", "operator_norm"):
        plus = np.asarray(observations[1][name])
        minus = np.asarray(observations[-1][name])
        paired[f"plus_minus_{name}"] = summary(plus - minus)

    return {
        "child_order": r,
        "parent_order": n,
        "source": source,
        "conference_identity_verified": True,
        "bridge_samples": samples,
        "independent_statistical_units": samples,
        "orientation_values_per_bridge": 2,
        "seed": seed,
        "theory": theoretical_moments(order),
        "orientations": output_orientations,
        "paired_orientation_differences": paired,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orders", type=int, nargs="+", default=sorted(SOURCES))
    parser.add_argument("--seed", type=int, default=20260814)
    parser.add_argument("--sample-multiplier", type=float, default=1.0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    unknown = sorted(set(args.orders) - set(SOURCES))
    if unknown:
        raise ValueError(f"no saved conference source for orders {unknown}")

    payload = {
        "schema": "conference-random-bridge-spectrum-v1",
        "classification": (
            "exact finite matrix powers for seeded Monte Carlo bridges; "
            "the expectation formulas through m6 are proved algebraically"
        ),
        "normalization": "W_eps=[[A,B],[B^T,eps*A]]/sqrt(2r)",
        "free_limit": (
            "Bernoulli(+-1/sqrt(2)) boxplus semicircle(variance=1/2)"
        ),
        "records": [
            audit_order(
                order,
                max(1, int(round(DEFAULT_SAMPLES[order] * args.sample_multiplier))),
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
