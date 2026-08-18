#!/usr/bin/env python3
"""Exact finite audit of the conference-parent quartic identity.

This is a task-local verifier for

    S = [[A, B], [B.T, epsilon A]],
    J = ||B B.T||_F^2 + ||A B + epsilon B A||_F^2.

It checks the trace/cumulant formula directly, and checks the spectral
projection representation of the second summand.  The order-two bridge
space is exhausted; order six uses a fixed seeded sample.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]


def conference(order: int) -> np.ndarray:
    if order == 2:
        return np.asarray([[0, 1], [1, 0]], dtype=np.int64)
    payload = json.loads(
        (ROOT / "computations" / "results" / f"exact_m{order}.json").read_text()
    )
    matrix = np.asarray(payload["matrix"], dtype=np.int64)
    if not np.array_equal(matrix @ matrix, (order - 1) * np.eye(order, dtype=np.int64)):
        raise AssertionError("input is not a symmetric conference signing")
    return matrix


def spins(order: int) -> np.ndarray:
    masks = np.arange(1 << order, dtype=np.uint64)[:, None]
    bits = (masks >> np.arange(order, dtype=np.uint64)) & 1
    return (1 - 2 * bits).astype(np.int64)


def bridge_records(order: int, samples: int, seed: int) -> list[dict]:
    a = conference(order)
    z = spins(2 * order)
    rng = np.random.default_rng(seed)
    if order == 2:
        masks = range(1 << (order * order))
    else:
        masks = range(samples)
    records: list[dict] = []
    powers = np.arange(order * order, dtype=np.uint64)
    for mask in masks:
        if order == 2:
            bits = ((np.uint64(mask) >> powers) & 1).astype(np.int64)
            b = (1 - 2 * bits).reshape(order, order)
        else:
            b = rng.choice(np.asarray([-1, 1], dtype=np.int64), (order, order))
        for epsilon in (-1, 1):
            s = np.block([[a, b], [b.T, epsilon * a]])
            f = int(np.sum((b @ b.T) ** 2))
            g = int(np.sum((a @ b + epsilon * b @ a) ** 2))
            j = f + g
            trace4 = int(np.trace(np.linalg.matrix_power(s, 4)))
            predicted_trace4 = (
                6 * order**3 - 8 * order**2 + 2 * order + 2 * j
            )
            if trace4 != predicted_trace4:
                raise AssertionError("trace-four identity failed")

            energies = np.einsum("bi,ij,bj->b", z, s, z, dtype=np.int64) // 2
            second = int(np.mean(energies**2))
            fourth = int(np.mean(energies**4))
            cumulant4 = fourth - 3 * second**2
            predicted_cumulant4 = (
                6 * j - 30 * order**3 + 32 * order**2 - 10 * order
            )
            if cumulant4 != predicted_cumulant4:
                raise AssertionError("fourth-cumulant identity failed")

            eigenvalues, u = np.linalg.eigh(a.astype(np.float64))
            transformed = u.T @ b @ u
            if epsilon == 1:
                active = np.abs(eigenvalues[:, None] + eigenvalues[None, :]) > 1e-8
            else:
                active = np.abs(eigenvalues[:, None] - eigenvalues[None, :]) > 1e-8
            projected_norm2 = float(np.sum(transformed[active] ** 2))
            predicted_g = 4 * (order - 1) * projected_norm2
            if not np.isclose(g, predicted_g, atol=1e-7):
                raise AssertionError("projection representation failed")
            records.append(
                {
                    "epsilon": epsilon,
                    "frobenius_gram": f,
                    "frobenius_intertwiner": g,
                    "J": j,
                    "trace_S4": trace4,
                    "kappa4": cumulant4,
                }
            )
    return records


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples-6", type=int, default=64)
    parser.add_argument("--seed", type=int, default=20260817)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    groups = []
    for order, samples in ((2, 16), (6, args.samples_6)):
        records = bridge_records(order, samples, args.seed + order)
        groups.append(
            {
                "order": order,
                "bridge_draws": len(records) // 2,
                "orientations_per_bridge": 2,
                "min_J": min(row["J"] for row in records),
                "max_J": max(row["J"] for row in records),
                "records_sha256": hashlib.sha256(
                    json.dumps(records, sort_keys=True).encode("utf-8")
                ).hexdigest(),
            }
        )
    payload = {
        "schema": "conference-quartic-identity-audit-v1",
        "classification": (
            "exact arithmetic trace and Boolean-moment checks; exhaustive at r=2; "
            "seeded identity regression at r=6"
        ),
        "seed": args.seed,
        "groups": groups,
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered)


if __name__ == "__main__":
    main()
