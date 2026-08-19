#!/usr/bin/env python3
"""Canonical-disorder restriction audit in the exact switching quotient.

Switching invariance reduces the order-n disorder cube from
``2^binom(n,2)`` signings to ``2^binom(n-1,2)`` root-gauged signings.  The
quotient makes the complete canonical law feasible through order eight and
preserves all entropy, cap, marginal-KL, and partition-defect terms after
adding the explicit gauge entropy.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np


def internal_edges(n: int) -> list[tuple[int, int]]:
    return list(itertools.combinations(range(1, n), 2))


def root_caps(n: int, chunk: int = 8192) -> np.ndarray:
    edges = internal_edges(n)
    d = len(edges)
    spin_masks = np.arange(1 << (n - 1), dtype=np.uint64)
    spins = np.ones((len(spin_masks), n), dtype=np.int8)
    for j in range(1, n):
        spins[:, j] = 1 - 2 * ((spin_masks >> (j - 1)) & 1).astype(np.int8)
    root_energy = np.sum(spins[:, 1:], axis=1, dtype=np.int16)
    if not edges:
        return np.asarray([np.max(np.abs(root_energy))], dtype=np.int16)
    products = np.stack(
        [spins[:, i] * spins[:, j] for i, j in edges], axis=1
    ).astype(np.int16)
    output = np.empty(1 << d, dtype=np.int16)
    bit_positions = np.arange(d, dtype=np.uint64)
    for start in range(0, 1 << d, chunk):
        stop = min(start + chunk, 1 << d)
        masks = np.arange(start, stop, dtype=np.uint64)
        signs = 1 - 2 * ((masks[:, None] >> bit_positions) & 1).astype(np.int16)
        energies = signs @ products.T + root_energy[None, :]
        output[start:stop] = np.max(np.abs(energies), axis=1)
    return output


def entropy(probability: np.ndarray) -> float:
    positive = probability[probability > 0]
    return float(-np.dot(positive, np.log(positive)))


def canonical_root(
    cap: np.ndarray, n: int, beta: float
) -> tuple[np.ndarray, float]:
    log_weight = -beta * math.sqrt(n) * cap.astype(float)
    shift = float(np.max(log_weight))
    weight = np.exp(log_weight - shift)
    root_log_z = shift + math.log(float(np.sum(weight)))
    return weight / float(np.sum(weight)), root_log_z + (n - 1) * math.log(2)


def projected_masks(n: int, m: int) -> np.ndarray:
    source_edges = internal_edges(n)
    target_position = {edge: bit for bit, edge in enumerate(internal_edges(m))}
    masks = np.arange(1 << len(source_edges), dtype=np.uint64)
    projected = np.zeros(len(masks), dtype=np.uint64)
    for source_bit, edge in enumerate(source_edges):
        if edge in target_position:
            projected |= ((masks >> source_bit) & 1) << target_position[edge]
    return projected.astype(np.int64)


def audit_pair(
    cap_by_order: dict[int, np.ndarray], n: int, m: int, beta: float
) -> dict[str, float | int]:
    cap_n, cap_m = cap_by_order[n], cap_by_order[m]
    mu_n, log_z_n = canonical_root(cap_n, n, beta)
    mu_m, log_z_m = canonical_root(cap_m, m, beta)
    marginal = np.bincount(
        projected_masks(n, m), weights=mu_n, minlength=len(cap_m)
    )
    l_n, l_m = n * (n - 1) // 2, m * (m - 1) // 2
    q = l_m / l_n
    h_n = entropy(mu_n) + (n - 1) * math.log(2)
    h_s = entropy(marginal) + (m - 1) * math.log(2)
    shearer = h_s / q - h_n
    e_n = float(np.dot(mu_n, cap_n))
    e_s = float(np.dot(marginal, cap_m))
    energy_excess = math.sqrt(m) * e_s / q - math.sqrt(n) * e_n
    positive = marginal > 0
    kl = float(
        np.dot(marginal[positive], np.log(marginal[positive] / mu_m[positive]))
    )
    residual = beta * energy_excess - shearer
    partition_defect = log_z_n - log_z_m / q
    return {
        "N": n,
        "m": m,
        "beta": beta,
        "shearer_slack": shearer,
        "energy_excess": energy_excess,
        "residual": residual,
        "marginal_kl": kl,
        "scaled_marginal_kl": kl / q,
        "partition_defect": partition_defect,
        "identity_error": residual - kl / q - partition_defect,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=8)
    parser.add_argument("--min-parent", type=int, default=7)
    parser.add_argument("--betas", type=float, nargs="+", default=[0.5, 1.0, 2.0])
    parser.add_argument("--chunk", type=int, default=8192)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    cap_by_order = {
        n: root_caps(n, args.chunk) for n in range(2, args.max_n + 1)
    }
    records = [
        audit_pair(cap_by_order, n, m, beta)
        for beta in args.betas
        for n in range(args.min_parent, args.max_n + 1)
        for m in range(max(2, math.ceil(n / 3)), math.floor(2 * n / 3) + 1)
    ]
    payload = {
        "max_n": args.max_n,
        "betas": args.betas,
        "records": records,
        "root_cap_histograms": {
            str(n): {
                str(int(value)): int(count)
                for value, count in zip(*np.unique(cap, return_counts=True))
            }
            for n, cap in cap_by_order.items()
        },
        "max_abs_identity_error": max(
            abs(float(row["identity_error"])) for row in records
        ),
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(rendered + "\n", encoding="utf-8")
    else:
        print(rendered)


if __name__ == "__main__":
    main()
