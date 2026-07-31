#!/usr/bin/env python3
"""Test adaptive Gibbs-gradient bridge construction on exact finite children.

At a partial bridge with energy E(x,y), adding edge (i,j) with sign c changes
the absolute partition function by

    Z_new / Z_old = cosh(gamma) * (1 + c r_ij tanh(gamma)),

where r_ij = <tau x_i y_j>.  The greedy rule selects an unrevealed edge of
maximum |r_ij| and takes c=-sign(r_ij).  Enumeration is exact up to floating
evaluation of Gibbs weights; the final Boolean cap is evaluated in integers.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from bridge_block_cpsat import load_matrix, one_copy_energies
from exact_mn_milp import exact_profile, stable_matrix_hash


def construct(a: np.ndarray, b: np.ndarray, sign_b: int, gamma: float) -> dict[str, object]:
    x, ea = one_copy_energies(a)
    y_projective, eb_projective = one_copy_energies(sign_b * b)
    x = x.astype(np.int64)
    # x may be projectivized globally, but y must retain both relative block
    # orientations.  Omitting -y would miss half of the parent spin states.
    y = np.concatenate((y_projective, -y_projective), axis=0).astype(np.int64)
    eb = np.concatenate((eb_projective, eb_projective), axis=0)
    energy = ea[:, None].astype(np.int64) + eb[None, :].astype(np.int64)
    bridge = np.zeros((len(a), len(b)), dtype=np.int64)
    unrevealed = np.ones_like(bridge, dtype=bool)
    increments: list[dict[str, object]] = []

    for step in range(bridge.size):
        scaled = gamma * energy.astype(float)
        shift = float(np.max(np.abs(scaled)))
        # Both numerator and denominator are evaluated after the same shift.
        denominator_weights = np.exp(-shift) * np.cosh(scaled)
        numerator_weights = np.exp(-shift) * np.sinh(scaled)
        denominator = float(np.sum(denominator_weights))
        correlations = (x.T @ numerator_weights @ y) / denominator
        masked = np.where(unrevealed, np.abs(correlations), -1.0)
        i, j = map(int, np.unravel_index(np.argmax(masked), masked.shape))
        r = float(correlations[i, j])
        edge_sign = -1 if r > 0 else 1
        exact_increment = float(
            np.log(np.cosh(gamma))
            + np.log1p(edge_sign * r * np.tanh(gamma))
        )
        bridge[i, j] = edge_sign
        unrevealed[i, j] = False
        energy = energy + edge_sign * x[:, i, None] * y[None, :, j]
        increments.append(
            {
                "step": step,
                "edge": [i, j],
                "sign": edge_sign,
                "correlation": r,
                "absolute_correlation": abs(r),
                "log_partition_increment": exact_increment,
            }
        )

    final_cap = int(np.max(np.abs(energy)))
    parent = np.block([[a, bridge], [bridge.T, sign_b * b]])
    verified_profile = exact_profile(parent)
    if verified_profile["M"] != final_cap:
        raise AssertionError((verified_profile["M"], final_cap))
    child_caps = [int(np.max(np.abs(ea))), int(np.max(np.abs(eb)))]
    b_defect = final_cap ** (2 / 3) - sum(cap ** (2 / 3) for cap in child_caps)
    direct_log_partition = float(
        np.log(np.mean(np.cosh(gamma * energy.astype(float))))
    )
    base = ea[:, None].astype(np.int64) + eb[None, :].astype(np.int64)
    base_log_partition = float(
        np.log(np.mean(np.cosh(gamma * base.astype(float))))
    )
    increment_sum = float(sum(row["log_partition_increment"] for row in increments))
    if abs(direct_log_partition - base_log_partition - increment_sum) > 1e-8:
        raise AssertionError((direct_log_partition, base_log_partition, increment_sum))
    return {
        "orders": [len(a), len(b)],
        "sign_b": sign_b,
        "gamma": gamma,
        "bridge": bridge.tolist(),
        "parent_matrix_sha256": stable_matrix_hash(parent),
        "verified_parent_profile": verified_profile,
        "child_caps": child_caps,
        "final_cap": final_cap,
        "b_scale_composition_defect": b_defect,
        "zero_defect_energy_target": math.pow(
            sum(cap ** (2 / 3) for cap in child_caps), 3 / 2
        ),
        "base_log_partition": base_log_partition,
        "final_log_partition": direct_log_partition,
        "bridge_log_partition_cost": increment_sum,
        "sum_absolute_correlations": float(
            sum(row["absolute_correlation"] for row in increments)
        ),
        "max_absolute_correlation": float(
            max(row["absolute_correlation"] for row in increments)
        ),
        "increments": increments,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orders", type=int, nargs="+", default=[4, 5, 6, 7, 8])
    parser.add_argument("--scaled-temperatures", type=float, nargs="+", default=[1, 2, 4, 8])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    cases = []
    for order in args.orders:
        child = load_matrix(Path(f"computations/results/exact_m{order}.json"))
        for scaled_temperature in args.scaled_temperatures:
            gamma = scaled_temperature / np.sqrt(2 * order)
            for sign_b in (1, -1):
                row = construct(child, child, sign_b, gamma)
                row["scaled_temperature"] = scaled_temperature
                cases.append(row)
                print(
                    f"{order}+{order} sign={sign_b:+d} t={scaled_temperature:g} "
                    f"cap={row['final_cap']} cost={row['bridge_log_partition_cost']:.6f} "
                    f"sum|r|={row['sum_absolute_correlations']:.6f}"
                )
    payload = {
        "schema": "quadratic-signing-greedy-gibbs-bridge-v1",
        "classification": "exact finite enumeration with floating Gibbs weights; heuristic scaling evidence",
        "rule": "reveal the unrevealed edge with largest absolute Gibbs correlation and choose the opposite sign",
        "cases": cases,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
