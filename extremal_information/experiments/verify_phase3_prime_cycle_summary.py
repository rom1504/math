#!/usr/bin/env python3
"""Exact finite checks for the prime-cycle response/congruence separation.

For several odd primes and distortions, this script constructs the mesh net
from Theorem CSC.2, computes its *actual* all-context sup error exhaustively,
and compares the distortion with the best one-state closed decoder.  The
classification of congruences of a prime cyclic group is proved separately;
the script only checks the finite response geometry and constants.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def response_distance(p: int, x: int, y: int) -> float:
    return max(
        abs(
            math.cos(2 * math.pi * (x + c) / p)
            - math.cos(2 * math.pi * (y + c) / p)
        )
        for c in range(p)
    )


def mesh_centers(p: int, epsilon: float) -> list[int]:
    k = math.ceil(2 * math.pi / epsilon)
    if p <= k:
        return list(range(p))
    return sorted({int(math.floor(j * p / k + 0.5)) % p for j in range(k)})


def run() -> dict[str, object]:
    cases: list[dict[str, object]] = []
    for p in (5, 11, 31, 101):
        values = [math.cos(2 * math.pi * x / p) for x in range(p)]
        one_state_error = (max(values) - min(values)) / 2
        formula_error = (1 + math.cos(math.pi / p)) / 2
        assert abs(one_state_error - formula_error) < 1e-12

        for epsilon in (0.2, 0.5, 0.9):
            centers = mesh_centers(p, epsilon)
            worst_net_error = max(
                min(response_distance(p, x, y) for y in centers)
                for x in range(p)
            )
            assert len(centers) <= math.ceil(2 * math.pi / epsilon) + 1
            assert worst_net_error <= epsilon + 1e-12
            cases.append(
                {
                    "p": p,
                    "epsilon": epsilon,
                    "center_count": len(centers),
                    "actual_uniform_net_error": worst_net_error,
                    "best_one_state_closed_error": one_state_error,
                    "closed_state_count": 1 if epsilon >= one_state_error else p,
                }
            )
    return {"status": "passed", "cases": cases}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name(
            "phase3_prime_cycle_summary_results.json"
        ),
    )
    args = parser.parse_args()
    result = run()
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
