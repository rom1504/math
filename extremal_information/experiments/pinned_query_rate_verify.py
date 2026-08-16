#!/usr/bin/env python3
"""Exhaustively verify the finite pinned-query decoder through order five."""

from __future__ import annotations

import argparse
import json
from itertools import combinations, product
from pathlib import Path


def verify_order(n: int) -> dict:
    edges = list(combinations(range(n), 2))
    edge_count = len(edges)
    energy_quantum = 1
    field_strength = n  # Strictly greater than a(n-1).
    spins = list(product((-1, 1), repeat=n))
    seen = set()

    for coefficients in product((-1, 1), repeat=edge_count):
        def quadratic(x: tuple[int, ...]) -> int:
            return sum(
                coefficients[k] * x[i] * x[j]
                for k, (i, j) in enumerate(edges)
            )

        maximum = max(quadratic(x) for x in spins)
        response = {
            u: max(
                quadratic(x)
                - maximum
                + field_strength * sum(ui * xi for ui, xi in zip(u, x))
                for x in spins
            )
            for u in spins
        }
        assert all(
            response[u]
            == field_strength * n + quadratic(u) - maximum
            for u in spins
        )

        recovered = []
        for i, j in edges:
            coefficient = sum(
                (response[u] - field_strength * n) * u[i] * u[j]
                for u in spins
            ) / (2**n)
            recovered.append(1 if coefficient > 0 else -1)
        assert tuple(recovered) == coefficients
        seen.add(tuple(response[u] for u in spins))

    assert len(seen) == 2**edge_count
    return {
        "order": n,
        "edges": edge_count,
        "landscapes_checked": 2**edge_count,
        "distinct_response_vectors": len(seen),
        "field_strength": field_strength,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-order", type=int, default=5)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "extremal_information/experiments/pinned_query_rate_results.json"
        ),
    )
    args = parser.parse_args()
    result = {
        "schema": "extremal-information-pinned-query-rate-v1",
        "orders": [verify_order(n) for n in range(2, args.max_order + 1)],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
