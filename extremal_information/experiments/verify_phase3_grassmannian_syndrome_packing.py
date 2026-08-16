#!/usr/bin/env python3
"""Exhaustive small checks for the Grassmannian syndrome-response packing.

The default check enumerates all subspaces of F_2^(2d) for d <= 2, selects a
canonical complement basis for each, and checks:

  * r(S_W) = d + 1;
  * r(S_W union S_W') <= dim(W intersection W') + 2;
  * the self-query response gap is at least d-dim(intersection)-1.

It uses integers as bit vectors and never writes outside the repository.
"""

from __future__ import annotations

import argparse
import json
from itertools import combinations
from pathlib import Path


def rank(vectors: tuple[int, ...] | list[int]) -> int:
    pivots: dict[int, int] = {}
    for value in vectors:
        x = value
        while x:
            bit = x.bit_length() - 1
            if bit in pivots:
                x ^= pivots[bit]
            else:
                pivots[bit] = x
                break
    return len(pivots)


def span(basis: tuple[int, ...] | list[int]) -> frozenset[int]:
    values = {0}
    for vector in basis:
        values |= {x ^ vector for x in tuple(values)}
    return frozenset(values)


def all_subspaces(w: int, d: int) -> list[frozenset[int]]:
    seen: set[frozenset[int]] = set()
    nonzero = range(1, 1 << w)
    for candidate in combinations(nonzero, d):
        if rank(candidate) == d:
            seen.add(span(candidate))
    return sorted(seen, key=lambda subspace: tuple(sorted(subspace)))


def canonical_complement_basis(W: frozenset[int], w: int) -> tuple[int, ...]:
    current = list(sorted(W))
    current_rank = rank(current)
    complement: list[int] = []
    for vector in range(1, 1 << w):
        new_rank = rank(current + [vector])
        if new_rank > current_rank:
            complement.append(vector)
            current.append(vector)
            current_rank = new_rank
            if current_rank == w:
                return tuple(complement)
    raise AssertionError("failed to find complement")


def cayley_radius(support: frozenset[int], w: int) -> int:
    distance = [w + 1] * (1 << w)
    distance[0] = 0
    frontier = [0]
    for depth in range(w + 1):
        next_frontier: list[int] = []
        for x in frontier:
            assert distance[x] == depth
            for generator in support:
                y = x ^ generator
                if distance[y] > depth + 1:
                    distance[y] = depth + 1
                    next_frontier.append(y)
        if not next_frontier:
            break
        frontier = next_frontier
    assert max(distance) <= w
    return max(distance)


def run(max_d: int) -> dict[str, object]:
    cases: list[dict[str, int]] = []
    for d in range(1, max_d + 1):
        w = 2 * d
        subspaces = all_subspaces(w, d)
        supports: dict[frozenset[int], frozenset[int]] = {}
        for W in subspaces:
            complement = canonical_complement_basis(W, w)
            support = frozenset((W - {0}) | set(complement))
            supports[W] = support
            assert cayley_radius(support, w) == d + 1

        pair_count = 0
        minimum_slack = w
        for i, W in enumerate(subspaces):
            for Wp in subspaces[i + 1 :]:
                t = (len(W & Wp)).bit_length() - 1
                assert len(W & Wp) == 1 << t
                cross = cayley_radius(supports[W] | supports[Wp], w)
                minimum_slack = min(minimum_slack, t + 2 - cross)
                assert cross <= t + 2
                self_gap = (d + 1) - cross
                assert self_gap >= d - t - 1
                pair_count += 1
        cases.append(
            {
                "d": d,
                "w": w,
                "subspace_count": len(subspaces),
                "pair_count": pair_count,
                "minimum_cross_bound_slack": minimum_slack,
            }
        )
    return {"status": "passed", "cases": cases}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-d", type=int, default=2)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name(
            "phase3_grassmannian_syndrome_packing_results.json"
        ),
    )
    args = parser.parse_args()
    result = run(args.max_d)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
