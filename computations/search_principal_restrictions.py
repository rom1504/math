#!/usr/bin/env python3
"""Exhaustively find the lowest-cap principal restrictions of fixed orders."""

from __future__ import annotations

import argparse
import itertools
import json
from collections import Counter
from pathlib import Path

import numpy as np

from exact_mn_milp import exact_profile, projective_spins, stable_matrix_hash


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("--matrix-key", required=True)
    parser.add_argument("--order", action="append", type=int, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = json.loads(args.source.read_text())
    matrix = np.asarray(payload[args.matrix_key], dtype=np.int8)
    records = []
    for order in sorted(set(args.order)):
        if not 2 <= order <= len(matrix):
            raise ValueError(order)
        spins = projective_spins(order).astype(np.int16)
        best_cap = 10**9
        best_subset = None
        best_matrix = None
        histogram: Counter[int] = Counter()
        checked = 0
        for subset in itertools.combinations(range(len(matrix)), order):
            principal = matrix[np.ix_(subset, subset)].astype(np.int16)
            energies = np.einsum("bi,ij,bj->b", spins, principal, spins) // 2
            cap = int(np.max(np.abs(energies)))
            histogram[cap] += 1
            checked += 1
            if cap < best_cap:
                best_cap = cap
                best_subset = subset
                best_matrix = principal.astype(np.int8)
        assert best_subset is not None and best_matrix is not None
        profile = exact_profile(best_matrix)
        if profile["M"] != best_cap:
            raise AssertionError((profile["M"], best_cap))
        records.append(
            {
                "order": order,
                "subsets_checked": checked,
                "cap_histogram": {
                    str(cap): count for cap, count in sorted(histogram.items())
                },
                "best_cap": best_cap,
                "best_subset": list(best_subset),
                "best_matrix": [
                    [int(value) for value in row] for row in best_matrix
                ],
                "best_matrix_sha256": stable_matrix_hash(best_matrix),
                "best_profile": profile,
            }
        )
        print(
            f"order={order} checked={checked} best_cap={best_cap} "
            f"histogram={dict(sorted(histogram.items()))}"
        )
    output = {
        "schema": "quadratic-signing-principal-restriction-search-v1",
        "classification": "exhaustive exact finite restriction search",
        "source": str(args.source),
        "matrix_key": args.matrix_key,
        "source_order": len(matrix),
        "source_matrix_sha256": stable_matrix_hash(matrix),
        "records": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
