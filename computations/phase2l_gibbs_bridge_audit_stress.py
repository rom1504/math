#!/usr/bin/env python3
"""Independent stable audit and stress test of the greedy Gibbs bridge.

Uses log-sum-exp normalized cosh weights, audits saved greedy trajectories,
and compares greedy with reproducible random reveal orders and held-out exact
minimizer classes.  All final caps are exact integer evaluations; Gibbs
quantities are floating evaluations of finite sums.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from bridge_block_cpsat import one_copy_energies


def log_cosh(values: np.ndarray) -> np.ndarray:
    absolute = np.abs(values)
    return absolute + np.log1p(np.exp(-2 * absolute)) - math.log(2.0)


def response(
    energy: np.ndarray, x: np.ndarray, y: np.ndarray, gamma: float
) -> tuple[float, np.ndarray]:
    scaled = gamma * energy.astype(float)
    log_weights = log_cosh(scaled)
    shift = float(np.max(log_weights))
    weights = np.exp(log_weights - shift)
    denominator = float(np.sum(weights))
    correlations = (x.T @ (weights * np.tanh(scaled)) @ y) / denominator
    log_partition = shift + math.log(denominator) - math.log(energy.size)
    return log_partition, correlations


def prepare(
    a: np.ndarray, b: np.ndarray, sign_b: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    x, ea = one_copy_energies(a)
    yp, ebp = one_copy_energies(sign_b * b)
    x = x.astype(np.int64)
    y = np.concatenate((yp, -yp), axis=0).astype(np.int64)
    eb = np.concatenate((ebp, ebp), axis=0).astype(np.int64)
    base = ea[:, None].astype(np.int64) + eb[None, :]
    return x, y, ea.astype(np.int64), eb, base


def run_strategy(
    a: np.ndarray,
    b: np.ndarray,
    sign_b: int,
    gamma: float,
    strategy: str,
    seed: int,
) -> dict[str, object]:
    x, y, ea, eb, base = prepare(a, b, sign_b)
    energy = base.copy()
    m, n = len(a), len(b)
    bridge = np.zeros((m, n), dtype=np.int64)
    unrevealed = np.ones((m, n), dtype=bool)
    rng = np.random.default_rng(seed)
    random_order = list(rng.permutation(m * n))
    absolute_correlations = []
    increment_sum = 0.0
    initial_log_partition, _ = response(energy, x, y, gamma)
    for step in range(m * n):
        current_log_partition, correlations = response(energy, x, y, gamma)
        if strategy == "greedy":
            masked = np.where(unrevealed, np.abs(correlations), -1.0)
            i, j = map(int, np.unravel_index(np.argmax(masked), masked.shape))
        elif strategy == "random":
            flat = int(random_order[step])
            i, j = divmod(flat, n)
            if not unrevealed[i, j]:
                raise AssertionError("random order repeated an edge")
        else:
            raise ValueError(strategy)
        r = float(correlations[i, j])
        edge_sign = -1 if r > 0 else 1
        increment = math.log(math.cosh(gamma)) + math.log1p(
            edge_sign * r * math.tanh(gamma)
        )
        bridge[i, j] = edge_sign
        unrevealed[i, j] = False
        energy = energy + edge_sign * x[:, i, None] * y[None, :, j]
        next_log_partition, _ = response(energy, x, y, gamma)
        if abs(next_log_partition - current_log_partition - increment) > 2e-12:
            raise AssertionError("increment identity failed")
        absolute_correlations.append(abs(r))
        increment_sum += increment
    final_log_partition, _ = response(energy, x, y, gamma)
    if abs(final_log_partition - initial_log_partition - increment_sum) > 2e-11:
        raise AssertionError("telescoping partition identity failed")
    child_caps = [int(np.max(np.abs(ea))), int(np.max(np.abs(eb)))]
    target = math.pow(sum(cap ** (2 / 3) for cap in child_caps), 3 / 2)
    delta = target - sum(child_caps)
    edges = m * n
    exact_linear_sufficient_threshold = (
        edges * math.log(math.cosh(gamma)) - gamma * delta
    ) / math.tanh(gamma)
    values = np.asarray(absolute_correlations)
    return {
        "orders": [m, n],
        "sign_b": sign_b,
        "gamma": gamma,
        "strategy": strategy,
        "seed": seed,
        "child_caps": child_caps,
        "zero_defect_energy_target": target,
        "delta": delta,
        "sum_absolute_correlations": float(np.sum(values)),
        "sum_absolute_correlations_over_n_three_halves": (
            float(np.sum(values)) / (m * n) ** 0.75
        ),
        "median_absolute_correlation": float(np.median(values)),
        "max_absolute_correlation": float(np.max(values)),
        "zero_correlation_step_count": int(np.count_nonzero(values < 1e-12)),
        "small_gamma_threshold": edges * gamma / 2 - delta,
        "exact_linear_sufficient_threshold": exact_linear_sufficient_threshold,
        "linear_threshold_surplus": (
            float(np.sum(values)) - exact_linear_sufficient_threshold
        ),
        "bridge_log_partition_cost": increment_sum,
        "calibrated_log_cost_margin": gamma * delta - increment_sum,
        "final_cap": int(np.max(np.abs(energy))),
        "bridge": bridge.tolist(),
    }


def audit_saved_case(
    row: dict[str, object], a: np.ndarray, b: np.ndarray
) -> dict[str, object]:
    x, y, _ea, _eb, energy = prepare(a, b, int(row["sign_b"]))
    gamma = float(row["gamma"])
    m, n = len(a), len(b)
    unrevealed = np.ones((m, n), dtype=bool)
    max_correlation_error = 0.0
    max_increment_error = 0.0
    max_greedy_choice_error = 0.0
    for saved in row["increments"]:
        before, correlations = response(energy, x, y, gamma)
        i, j = map(int, saved["edge"])
        if not unrevealed[i, j]:
            raise AssertionError("saved trajectory repeated an edge")
        masked = np.where(unrevealed, np.abs(correlations), -1.0)
        max_available = float(np.max(masked))
        r = float(correlations[i, j])
        sign = int(saved["sign"])
        expected_sign = -1 if r > 0 else 1
        # At an exact zero gradient either sign is a minimizer.  Independent
        # floating summation can flip signs at the 1e-17 scale.
        if abs(r) > 1e-12 and sign != expected_sign:
            raise AssertionError("saved trajectory did not oppose the Gibbs gradient")
        max_greedy_choice_error = max(
            max_greedy_choice_error, max_available - abs(r)
        )
        unrevealed[i, j] = False
        energy = energy + sign * x[:, i, None] * y[None, :, j]
        after, _ = response(energy, x, y, gamma)
        max_correlation_error = max(
            max_correlation_error, abs(r - float(saved["correlation"]))
        )
        max_increment_error = max(
            max_increment_error,
            abs(after - before - float(saved["log_partition_increment"])),
        )
    if (
        max_correlation_error > 2e-12
        or max_increment_error > 2e-12
        or max_greedy_choice_error > 2e-12
    ):
        raise AssertionError(
            (
                max_correlation_error,
                max_increment_error,
                max_greedy_choice_error,
            )
        )
    if int(np.max(np.abs(energy))) != int(row["final_cap"]):
        raise AssertionError("saved final cap changed")
    return {
        "orders": row["orders"],
        "sign_b": row["sign_b"],
        "gamma": gamma,
        "max_correlation_error": max_correlation_error,
        "max_increment_error": max_increment_error,
        "max_greedy_choice_error": max_greedy_choice_error,
        "final_cap": row["final_cap"],
        "verified": True,
    }


def load_exact(order: int) -> np.ndarray:
    payload = json.loads(
        Path(f"computations/results/exact_m{order}.json").read_text()
    )
    return np.asarray(payload["matrix"], dtype=np.int64)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--random-seeds", type=int, default=5)
    args = parser.parse_args()

    saved_rows = []
    for path in (
        Path("computations/results/phase2k_greedy_gibbs_bridge.json"),
        Path("computations/results/phase2k_greedy_gibbs_bridge_m9.json"),
        Path("computations/results/phase2k_greedy_gibbs_bridge_m10.json"),
    ):
        saved_rows.extend(json.loads(path.read_text())["cases"])
    audits = []
    for order in range(4, 11):
        candidates = [
            row for row in saved_rows
            if row["orders"] == [order, order]
            and row["sign_b"] == 1
            and row.get("scaled_temperature") == 4
        ]
        if len(candidates) != 1:
            raise AssertionError((order, len(candidates)))
        child = load_exact(order)
        audits.append(audit_saved_case(candidates[0], child, child))

    stress = []
    for order in range(4, 9):
        child = load_exact(order)
        for scaled_temperature in (2.0, 4.0):
            gamma = scaled_temperature / math.sqrt(2 * order)
            stress.append(run_strategy(child, child, 1, gamma, "greedy", 0))
            for seed in range(args.random_seeds):
                stress.append(
                    run_strategy(child, child, 1, gamma, "random", seed)
                )

    heldout_classes = []
    for order in (7, 8):
        payload = json.loads(
            Path(f"computations/results/m{order}_minimizer_orbits.json").read_text()
        )
        gamma = 4.0 / math.sqrt(2 * order)
        for row in payload["classes"]:
            child = np.asarray(row["representative_matrix"], dtype=np.int64)
            result = run_strategy(child, child, 1, gamma, "greedy", 0)
            result["child_class"] = row["class"]
            result["child_class_orbit_sha256"] = row["canonical_orbit_sha256"]
            heldout_classes.append(result)

    payload = {
        "schema": "quadratic-signing-gibbs-bridge-stable-audit-v1",
        "classification": (
            "independent stable finite-sum audit and reproducible strategy stress; "
            "heuristic scaling evidence, no uniform gradient theorem"
        ),
        "normalization": (
            "x is projective and y includes both orientations, giving exactly "
            "2^(m+n-1) parent states; weights are proportional to cosh(gamma E)"
        ),
        "saved_trajectory_audits": audits,
        "strategy_stress": stress,
        "heldout_minimizer_classes": heldout_classes,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    print(f"audited {len(audits)} saved trajectories")
    print(f"strategy stress cases={len(stress)} heldout classes={len(heldout_classes)}")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
