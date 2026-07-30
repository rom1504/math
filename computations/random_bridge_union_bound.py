#!/usr/bin/env python3
"""Rigorous finite random-bridge union bound for two fixed child signings.

For fixed child Hamiltonians and an iid random sign bridge C, each x^T C y
has the law of a sum S_k of k=m*n independent signs.  After accounting for
the relative global spin flip between the blocks, the assembled cap is

    max_(x,y projective) (|H_A(x) + s H_B(y)| + |x^T C y|).

Thus a bridge of cap at most T exists whenever

    sum_(x,y) Pr(|S_(mn)| + |H_A(x)+s H_B(y)| > T) < 1.

The script evaluates the numerator of this union bound exactly as an integer
over the common denominator 2^(mn), using energy histograms.  This is a
rigorous finite existence certificate, although it need not be sharp.
"""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np

from exact_mn_milp import projective_spins


def load_matrix(path: Path) -> np.ndarray:
    payload = json.loads(path.read_text())
    key = "matrix" if "matrix" in payload else "parent_matrix"
    matrix = np.asarray(payload[key], dtype=np.int8)
    if not np.array_equal(matrix, matrix.T) or np.any(np.diag(matrix)):
        raise ValueError(f"invalid signing matrix in {path}")
    return matrix


def energy_histogram(matrix: np.ndarray) -> Counter[int]:
    spins = projective_spins(len(matrix)).astype(np.int64)
    energy = np.einsum("bi,ij,bj->b", spins, matrix.astype(np.int64), spins) // 2
    return Counter(int(v) for v in energy)


def internal_absolute_histogram(
    hist_a: Counter[int], hist_b: Counter[int], sign_b: int
) -> Counter[int]:
    result: Counter[int] = Counter()
    for ea, count_a in hist_a.items():
        for eb, count_b in hist_b.items():
            result[abs(ea + sign_b * eb)] += count_a * count_b
    return result


def strict_absolute_tail_numerator(k: int, threshold: int) -> int:
    """Return 2^k Pr(|S_k|>threshold) exactly as an integer."""
    if threshold < 0:
        return 1 << k
    if threshold >= k:
        return 0
    numerator = 0
    for negative_count in range(k + 1):
        value = k - 2 * negative_count
        if abs(value) > threshold:
            numerator += math.comb(k, negative_count)
    return numerator


def union_numerator(k: int, internal_hist: Counter[int], cap: int) -> int:
    return sum(
        state_count * strict_absolute_tail_numerator(k, cap - internal)
        for internal, state_count in internal_hist.items()
    )


def allowed_ceiling(value: float, parity: int) -> int:
    cap = math.ceil(value - 1e-12)
    if cap % 2 != parity:
        cap += 1
    return cap


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("child_a", type=Path)
    parser.add_argument("child_b", type=Path)
    parser.add_argument("--sign-b", type=int, choices=(-1, 1), default=1)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    a = load_matrix(args.child_a)
    b = load_matrix(args.child_b)
    m, n = len(a), len(b)
    hist_a = energy_histogram(a)
    hist_b = energy_histogram(b)
    internal_hist = internal_absolute_histogram(hist_a, hist_b, args.sign_b)
    bridge_edges = m * n
    denominator = 1 << bridge_edges
    total_edges = (m + n) * (m + n - 1) // 2
    parity = total_edges % 2
    cap_a = max(abs(v) for v in hist_a)
    cap_b = max(abs(v) for v in hist_b)
    ideal = (cap_a ** (2.0 / 3.0) + cap_b ** (2.0 / 3.0)) ** 1.5
    ideal_cap = allowed_ceiling(ideal, parity)

    minimum_cap = None
    minimum_numerator = None
    audit_rows = []
    start = min(internal_hist)
    if start % 2 != parity:
        start += 1
    for cap in range(start, total_edges + 1, 2):
        numerator = union_numerator(bridge_edges, internal_hist, cap)
        if cap in {ideal_cap - 2, ideal_cap, ideal_cap + 2}:
            audit_rows.append(
                {
                    "cap": cap,
                    "union_numerator": str(numerator),
                    "union_bound": numerator / denominator,
                    "certifies_existence": numerator < denominator,
                }
            )
        if numerator < denominator:
            minimum_cap = cap
            minimum_numerator = numerator
            break

    payload = {
        "schema": "quadratic-signing-random-bridge-union-v1",
        "classification": "proved finite existence bound by exact integer union-bound arithmetic",
        "child_a": str(args.child_a),
        "child_b": str(args.child_b),
        "orders": [m, n],
        "sign_b": args.sign_b,
        "child_caps": [cap_a, cap_b],
        "bridge_random_variables": bridge_edges,
        "union_denominator": str(denominator),
        "internal_absolute_histogram": {
            str(value): count for value, count in sorted(internal_hist.items())
        },
        "ideal_two_thirds_energy_target": ideal,
        "ideal_allowed_cap": ideal_cap,
        "audit_rows": audit_rows,
        "minimum_certified_cap": minimum_cap,
        "minimum_union_numerator": (
            None if minimum_numerator is None else str(minimum_numerator)
        ),
        "minimum_union_bound": (
            None if minimum_numerator is None else minimum_numerator / denominator
        ),
    }
    print(
        f"{m}+{n} sign_b={args.sign_b:+d}: ideal={ideal:.12f} "
        f"allowed={ideal_cap} union-certified-cap={minimum_cap} "
        f"bound={payload['minimum_union_bound']}"
    )
    for row in audit_rows:
        print(
            f"  cap={row['cap']} union={row['union_bound']:.12g} "
            f"exists={row['certifies_existence']}"
        )
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
        print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
