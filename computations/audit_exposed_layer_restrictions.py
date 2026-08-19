#!/usr/bin/env python3
"""Exact restriction tails for finite-temperature exposed cap layers.

Switching quotient enumeration is exact through order eight.  For each
``(n, beta)`` this script selects a cumulative cap level maximizing

    L_n^{-1} log |{A: Q(A) <= c n^(3/2)}| - beta n^2 c / L_n

and computes the normalized-cap tail of a uniform proportional principal
restriction.  These are the finite quantities in (ERSR) of
``artifacts/cross_order_outward_director_review.md``.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np

from audit_canonical_disorder_root_gauge import projected_masks, root_caps


def exposed_level(cap: np.ndarray, n: int, beta: float) -> tuple[int, float, int]:
    values, counts = np.unique(cap, return_counts=True)
    edge_count = n * (n - 1) // 2
    coefficient = beta * n * n / edge_count
    cumulative = 0
    best: tuple[float, int, float, int] | None = None
    for value, count in zip(values, counts):
        cumulative += int(count)
        normalized = float(value) / n**1.5
        score = math.log(cumulative) / edge_count - coefficient * normalized
        candidate = (score, int(value), normalized, cumulative)
        if best is None or candidate[0] > best[0]:
            best = candidate
    assert best is not None
    return best[1], best[2], best[3]


def audit_pair(
    cap_by_order: dict[int, np.ndarray], n: int, m: int, beta: float
) -> dict[str, float | int]:
    parent_cap = cap_by_order[n]
    child_cap = cap_by_order[m]
    level, normalized_level, layer_size = exposed_level(parent_cap, n, beta)
    parent_law = (parent_cap <= level).astype(float)
    parent_law /= float(np.sum(parent_law))
    marginal = np.bincount(
        projected_masks(n, m), weights=parent_law, minlength=len(child_cap)
    )
    normalized_child = child_cap / m**1.5
    bad = normalized_child > normalized_level + 1e-12
    positive_gaps = normalized_child[bad] - normalized_level
    bad_mass = min(1.0, max(0.0, float(np.sum(marginal[bad]))))
    return {
        "N": n,
        "m": m,
        "beta": beta,
        "exposed_cap": level,
        "exposed_normalized_cap": normalized_level,
        "exposed_root_layer_size": layer_size,
        "bad_restriction_mass": bad_mass,
        "minimum_positive_normalized_gap": (
            float(np.min(positive_gaps)) if len(positive_gaps) else None
        ),
        "expected_restriction_normalized_cap": float(
            np.dot(marginal, normalized_child)
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-parent", type=int, default=7)
    parser.add_argument("--max-n", type=int, default=8)
    parser.add_argument("--betas", type=float, nargs="+", default=[0.25, 0.5, 1, 2])
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
    rendered = json.dumps({"records": records}, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(rendered + "\n", encoding="utf-8")
    else:
        print(rendered)


if __name__ == "__main__":
    main()
