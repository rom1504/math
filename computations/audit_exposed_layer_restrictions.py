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


def exposed_level(
    cap: np.ndarray, n: int, beta: float
) -> tuple[int, float, int, int, int | None]:
    values, counts = np.unique(cap, return_counts=True)
    edge_count = n * (n - 1) // 2
    coefficient = beta * n * n / edge_count
    cumulative = 0
    best: tuple[float, int, float, int, int, int | None] | None = None
    previous_value: int | None = None
    for value, count in zip(values, counts):
        previous_cumulative = cumulative
        cumulative += int(count)
        normalized = float(value) / n**1.5
        score = math.log(cumulative) / edge_count - coefficient * normalized
        candidate = (
            score,
            int(value),
            normalized,
            cumulative,
            previous_cumulative,
            previous_value,
        )
        if best is None or candidate[0] > best[0]:
            best = candidate
        previous_value = int(value)
    assert best is not None
    return best[1], best[2], best[3], best[4], best[5]


def audit_pair(
    cap_by_order: dict[int, np.ndarray], n: int, m: int, beta: float
) -> dict[str, float | int]:
    parent_cap = cap_by_order[n]
    child_cap = cap_by_order[m]
    level, normalized_level, layer_size, lower_size, lower_level = exposed_level(
        parent_cap, n, beta
    )
    parent_law = (parent_cap <= level).astype(float)
    parent_law /= float(np.sum(parent_law))
    lower_mass = lower_size / layer_size
    shell_bound = (
        math.exp(-beta * math.sqrt(n) * (level - lower_level))
        if lower_level is not None
        else 0.0
    )
    if lower_mass > shell_bound + 1e-12:
        raise AssertionError((n, beta, lower_mass, shell_bound))
    marginal = np.bincount(
        projected_masks(n, m), weights=parent_law, minlength=len(child_cap)
    )
    shell_law = (parent_cap == level).astype(float)
    shell_law /= float(np.sum(shell_law))
    shell_marginal = np.bincount(
        projected_masks(n, m), weights=shell_law, minlength=len(child_cap)
    )
    normalized_child = child_cap / m**1.5
    bad = normalized_child > normalized_level + 1e-12
    positive_gaps = normalized_child[bad] - normalized_level
    bad_mass = min(1.0, max(0.0, float(np.sum(marginal[bad]))))
    shell_bad_mass = min(1.0, max(0.0, float(np.sum(shell_marginal[bad]))))
    return {
        "N": n,
        "m": m,
        "beta": beta,
        "exposed_cap": level,
        "exposed_normalized_cap": normalized_level,
        "exposed_root_layer_size": layer_size,
        "lower_cumulative_root_size": lower_size,
        "lower_layer_mass": lower_mass,
        "top_shell_gap": level - lower_level if lower_level is not None else None,
        "top_shell_exponential_bound": shell_bound,
        "bad_restriction_mass": bad_mass,
        "exact_shell_bad_restriction_mass": shell_bad_mass,
        "minimum_positive_normalized_gap": (
            float(np.min(positive_gaps)) if len(positive_gaps) else None
        ),
        "expected_restriction_normalized_cap": float(
            np.dot(marginal, normalized_child)
        ),
        "exact_shell_expected_restriction_normalized_cap": float(
            np.dot(shell_marginal, normalized_child)
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
