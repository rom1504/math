#!/usr/bin/env python3
"""Reproducible finite witness search for the RD.1 margin mechanism.

This is diagnostic evidence only; the asymptotic theorem is probabilistic.
All arithmetic in the reported certificates is exact integer arithmetic.
"""

import argparse
import json
import random


def matvec(B, y):
    return [sum(a * b for a, b in zip(row, y)) for row in B]


def dot(x, z):
    return sum(a * b for a, b in zip(x, z))


def sign(z):
    return tuple(1 if v >= 0 else -1 for v in z)


def trial(n, queries, rng):
    B = [[rng.choice((-1, 1)) for _ in range(n)] for _ in range(n)]
    ys = [tuple(rng.choice((-1, 1)) for _ in range(n)) for _ in range(queries)]
    zs = [matvec(B, y) for y in ys]
    xs = [sign(z) for z in zs]
    diagonals = [dot(x, z) for x, z in zip(xs, zs)]
    gaps = []
    for c, z in enumerate(zs):
        competitor = max(dot(xs[d], z) for d in range(queries) if d != c)
        gaps.append(diagonals[c] - competitor)
    linear_separations = [
        sum(abs(a - b) for a, b in zip(zs[c], zs[d]))
        for c in range(queries)
        for d in range(c + 1, queries)
    ]
    return {
        "B": B,
        "queries": ys,
        "states": xs,
        "diagonals": diagonals,
        "minimum_gap": min(gaps),
        "minimum_linear_response_separation": min(linear_separations),
        "all_states_distinct": len(set(xs)) == queries,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=24)
    parser.add_argument("--queries", type=int, default=24)
    parser.add_argument("--trials", type=int, default=200)
    parser.add_argument("--seed", type=int, default=1729)
    args = parser.parse_args()
    rng = random.Random(args.seed)
    best = None
    for _ in range(args.trials):
        candidate = trial(args.n, args.queries, rng)
        if best is None or candidate["minimum_gap"] > best["minimum_gap"]:
            best = candidate
    certificate = {
        "n": args.n,
        "query_count": args.queries,
        "trials": args.trials,
        "seed": args.seed,
        "minimum_gap": best["minimum_gap"],
        "minimum_linear_response_separation":
            best["minimum_linear_response_separation"],
        "normalized_linear_separation_n_3_over_2":
            best["minimum_linear_response_separation"] / args.n ** 1.5,
        "normalized_gap_n_3_over_2": best["minimum_gap"] / args.n ** 1.5,
        "minimum_diagonal": min(best["diagonals"]),
        "all_states_distinct": best["all_states_distinct"],
        "B": best["B"],
        "queries": best["queries"],
        "states": best["states"],
        "diagonals": best["diagonals"],
    }
    print(json.dumps(certificate, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
