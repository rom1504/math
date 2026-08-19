#!/usr/bin/env python3
"""Exact small-order audit of canonical-disorder restriction compensation.

The program enumerates every hollow symmetric signing through ``--max-n``
(six is quick), computes its Boolean cap, and verifies Proposition DR.1 in
``artifacts/cross_order_outward_director_review.md``.  It also checks whether
the restriction of the uniform minimizing fibre is the uniform minimizing
fibre at the smaller order.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np


def edges(n: int) -> list[tuple[int, int]]:
    return [(i, j) for i in range(n) for j in range(i + 1, n)]


def caps(n: int, chunk: int = 1 << 15) -> np.ndarray:
    """Return Q(A) for all signings, indexed by the edge bit mask."""
    es = edges(n)
    e = len(es)
    spin_masks = np.arange(1 << (n - 1), dtype=np.uint64)
    spins = np.ones((len(spin_masks), n), dtype=np.int8)
    for j in range(1, n):
        spins[:, j] = 1 - 2 * ((spin_masks >> (j - 1)) & 1).astype(np.int8)
    characters = np.stack(
        [spins[:, i] * spins[:, j] for i, j in es], axis=1
    ).astype(np.int16)

    out = np.empty(1 << e, dtype=np.int16)
    for start in range(0, 1 << e, chunk):
        stop = min(start + chunk, 1 << e)
        masks = np.arange(start, stop, dtype=np.uint64)
        signs = np.empty((stop - start, e), dtype=np.int16)
        for j in range(e):
            signs[:, j] = 1 - 2 * ((masks >> j) & 1).astype(np.int16)
        energies = characters @ signs.T
        out[start:stop] = np.max(np.abs(energies), axis=0)
    return out


def entropy(p: np.ndarray) -> float:
    positive = p[p > 0]
    return float(-np.dot(positive, np.log(positive)))


def restriction_index_map(n: int, m: int) -> np.ndarray:
    """Map every order-n edge mask to its restriction on vertices [m]."""
    full_edges = edges(n)
    small_position = {edge: k for k, edge in enumerate(edges(m))}
    masks = np.arange(1 << len(full_edges), dtype=np.uint64)
    result = np.zeros(len(masks), dtype=np.uint64)
    for full_k, edge in enumerate(full_edges):
        if edge in small_position:
            bit = (masks >> full_k) & 1
            result |= bit << small_position[edge]
    return result.astype(np.int64)


def canonical(c: np.ndarray, n: int, beta: float) -> tuple[np.ndarray, float]:
    logw = -beta * math.sqrt(n) * c.astype(float)
    shift = float(np.max(logw))
    w = np.exp(logw - shift)
    z = float(np.sum(w))
    return w / z, shift + math.log(z)


def audit_pair(
    all_caps: dict[int, np.ndarray], n: int, m: int, beta: float
) -> dict[str, float | int | bool]:
    cap_n, cap_m = all_caps[n], all_caps[m]
    mu_n, logz_n = canonical(cap_n, n, beta)
    mu_m, logz_m = canonical(cap_m, m, beta)
    index = restriction_index_map(n, m)
    marginal = np.bincount(index, weights=mu_n, minlength=len(cap_m))

    l_n = n * (n - 1) // 2
    l_m = m * (m - 1) // 2
    q = l_m / l_n
    h_n, h_s = entropy(mu_n), entropy(marginal)
    shearer = h_s / q - h_n
    e_n = float(np.dot(mu_n, cap_n))
    e_s = float(np.dot(marginal, cap_m))
    energy_excess = math.sqrt(m) * e_s / q - math.sqrt(n) * e_n
    positive = marginal > 0
    kl = float(np.dot(marginal[positive], np.log(marginal[positive] / mu_m[positive])))
    residual = beta * energy_excess - shearer
    partition_defect = logz_n - logz_m / q
    identity_error = residual - kl / q - partition_defect

    min_n, min_m = int(np.min(cap_n)), int(np.min(cap_m))
    frozen_n = (cap_n == min_n).astype(float)
    frozen_n /= float(np.sum(frozen_n))
    frozen_m = (cap_m == min_m).astype(float)
    frozen_m /= float(np.sum(frozen_m))
    frozen_marginal = np.bincount(index, weights=frozen_n, minlength=len(cap_m))
    frozen_exact = bool(np.max(np.abs(frozen_marginal - frozen_m)) < 1e-14)

    return {
        "N": n,
        "m": m,
        "beta": beta,
        "shearer_slack": shearer,
        "energy_excess": energy_excess,
        "residual": residual,
        "marginal_kl": kl,
        "partition_defect": partition_defect,
        "identity_error": identity_error,
        "frozen_marginal_is_uniform_minimizer": frozen_exact,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-n", type=int, default=6)
    parser.add_argument("--betas", type=float, nargs="+", default=[0.5, 1.0, 2.0])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    all_caps = {n: caps(n) for n in range(2, args.max_n + 1)}
    records = [
        audit_pair(all_caps, n, m, beta)
        for beta in args.betas
        for n in range(3, args.max_n + 1)
        for m in range(2, n)
    ]
    payload = {
        "max_n": args.max_n,
        "betas": args.betas,
        "cap_histograms": {
            str(n): {
                str(int(value)): int(count)
                for value, count in zip(*np.unique(c, return_counts=True))
            }
            for n, c in all_caps.items()
        },
        "records": records,
        "max_abs_identity_error": max(abs(float(r["identity_error"])) for r in records),
    }
    rendered = json.dumps(payload, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(rendered + "\n", encoding="utf-8")
    else:
        print(rendered)


if __name__ == "__main__":
    main()
