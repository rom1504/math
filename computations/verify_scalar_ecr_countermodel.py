#!/usr/bin/env python3
"""Numerically audit the log-periodic scalar ECR countermodel.

The proof is analytic.  This script checks the parity rounding, eventual
one-vertex inequalities, and the fixed-ratio asymptotics in SC.8--SC.12 of
artifacts/cross_order_scalar_entropy_restriction_no_go.md.
"""

from __future__ import annotations

import argparse
import json
import math


LOG2 = math.log(2.0)
OMEGA = math.pi / LOG2


def edge_count(n: int) -> int:
    return n * (n - 1) // 2


def profile(x: float) -> float:
    return 0.40 + 0.01 * math.sin(OMEGA * math.log(x))


def raw_cap(n: int) -> float:
    return n**1.5 * profile(float(n))


def parity_round(value: float, parity: int) -> int:
    """Nearest integer to value congruent to parity modulo two."""
    return parity + 2 * math.floor((value - parity) / 2.0 + 0.5)


def cap(n: int) -> int:
    return parity_round(raw_cap(n), edge_count(n) & 1)


def phase_record(j: int, beta: float) -> dict[str, float | int]:
    n = round(4**j / math.sqrt(2.0))
    m = n // 2
    l_n = edge_count(n)
    l_m = edge_count(m)
    q = l_m / l_n
    d_excess = math.sqrt(m) * cap(m) / q - math.sqrt(n) * cap(n)
    psi_n = LOG2 - beta * math.sqrt(n) * cap(n) / l_n
    psi_m = LOG2 - beta * math.sqrt(m) * cap(m) / l_m
    phase = (OMEGA * math.log(n) + math.pi / 2.0) % (2.0 * math.pi)
    phase_error = min(phase, 2.0 * math.pi - phase)
    return {
        "j": j,
        "N": n,
        "m": m,
        "phase_error": phase_error,
        "T_N_over_N_3_2": cap(n) / n**1.5,
        "T_m_over_m_3_2": cap(m) / m**1.5,
        "D_over_N2": d_excess / n**2,
        "beta_D_over_N2": beta * d_excess / n**2,
        "psi_N_minus_psi_m": psi_n - psi_m,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-n", type=int, default=100)
    parser.add_argument("--max-n", type=int, default=100_000)
    parser.add_argument("--min-j", type=int, default=4)
    parser.add_argument("--max-j", type=int, default=10)
    parser.add_argument("--beta", type=float, default=1.0)
    args = parser.parse_args()

    increments = [cap(n + 1) - cap(n) for n in range(args.min_n, args.max_n)]
    violations = [
        n
        for n, increment in zip(range(args.min_n, args.max_n), increments)
        if not (0 < increment < n)
    ]
    report = {
        "scan": {
            "min_n": args.min_n,
            "max_n": args.max_n,
            "minimum_increment": min(increments),
            "maximum_increment": max(increments),
            "one_vertex_violations": violations[:20],
            "violation_count": len(violations),
        },
        "phase_subsequence": [
            phase_record(j, args.beta) for j in range(args.min_j, args.max_j + 1)
        ],
        "predicted_limits": {
            "D_over_N2": 0.02,
            "beta_D_over_N2": 0.02 * args.beta,
            "psi_N_minus_psi_m": 0.04 * args.beta,
        },
    }
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
