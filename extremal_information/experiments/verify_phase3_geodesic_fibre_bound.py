#!/usr/bin/env python3
"""Verify the geodesic-fibre structure for binary syndrome supports.

The exhaustive portion checks every spanning support in F_2^w for w <= 4,
every diametral vertex, and every shortest generator representation of that
vertex.  A deterministic random portion checks larger ranks.  The script
also records the sharp D=2 endpoint example where an anticode fibre can have
four, rather than D+1=3, elements.

No external solver or package is required.
"""

from __future__ import annotations

import argparse
import itertools
import json
import random
from collections import Counter, deque
from pathlib import Path


def rank_binary(vectors: tuple[int, ...] | list[int], w: int) -> int:
    pivots = [0] * w
    rank = 0
    for vector in vectors:
        x = vector
        while x:
            i = x.bit_length() - 1
            if pivots[i]:
                x ^= pivots[i]
            else:
                pivots[i] = x
                rank += 1
                break
    return rank


def word_lengths(support: tuple[int, ...], w: int) -> list[int]:
    size = 1 << w
    distance = [-1] * size
    distance[0] = 0
    queue = deque([0])
    while queue:
        x = queue.popleft()
        for s in support:
            y = x ^ s
            if distance[y] < 0:
                distance[y] = distance[x] + 1
                queue.append(y)
    return distance


def coordinates_in_basis(x: int, basis: tuple[int, ...]) -> int:
    d = len(basis)
    for mask in range(1 << d):
        total = 0
        for i, b in enumerate(basis):
            if mask >> i & 1:
                total ^= b
        if total == x:
            return mask
    raise AssertionError("vector is not in the declared span")


def span(basis: tuple[int, ...]) -> set[int]:
    values = {0}
    for b in basis:
        values |= {x ^ b for x in tuple(values)}
    return values


def verify_geodesic(
    support: tuple[int, ...], basis: tuple[int, ...], diameter: int
) -> tuple[bool, str]:
    if rank_binary(basis, max(support).bit_length()) != diameter:
        return False, "geodesic generators are dependent"
    w_space = span(basis)
    if set(support) & w_space != set(basis):
        return False, "zero fibre is not exactly the geodesic basis"

    coordinates = {x: coordinates_in_basis(x, basis) for x in w_space}
    fibres: dict[int, list[int]] = {}
    unused = set(support) - w_space
    while unused:
        first = min(unused)
        fibre = sorted(s for s in unused if (s ^ first) in w_space)
        for s in fibre:
            unused.remove(s)
        fibre_coordinates = [coordinates[s ^ first] for s in fibre]
        fibres[first] = fibre_coordinates

    for fibre in fibres.values():
        for x, y in itertools.combinations(fibre, 2):
            if bin(x ^ y).count("1") > 2:
                return False, "a fibre has Hamming diameter greater than two"
        if diameter >= 3 and len(fibre) > diameter + 1:
            return False, "a diameter-two fibre has more than D+1 points"
    return True, ""


def check_support(support: tuple[int, ...], w: int) -> tuple[int, int]:
    distance = word_lengths(support, w)
    if min(distance) < 0:
        raise AssertionError("support does not span")
    diameter = max(distance)
    checked = 0
    for target, target_distance in enumerate(distance):
        if target_distance != diameter:
            continue
        for basis in itertools.combinations(support, diameter):
            total = 0
            for b in basis:
                total ^= b
            if total != target:
                continue
            ok, reason = verify_geodesic(support, basis, diameter)
            if not ok:
                raise AssertionError(
                    f"w={w}, S={support}, t={target}, B={basis}: {reason}"
                )
            checked += 1
    if checked == 0:
        raise AssertionError("no shortest representation was found")
    return diameter, checked


def exhaustive(max_w: int) -> dict[str, object]:
    report: dict[str, object] = {}
    for w in range(1, max_w + 1):
        universe = tuple(range(1, 1 << w))
        support_count = 0
        geodesic_count = 0
        diameter_histogram: Counter[int] = Counter()
        for mask in range(1 << len(universe)):
            support = tuple(
                universe[i] for i in range(len(universe)) if mask >> i & 1
            )
            if rank_binary(support, w) != w:
                continue
            diameter, checked = check_support(support, w)
            support_count += 1
            geodesic_count += checked
            diameter_histogram[diameter] += 1
        report[str(w)] = {
            "spanning_supports": support_count,
            "diametral_geodesics_checked": geodesic_count,
            "diameter_histogram": dict(sorted(diameter_histogram.items())),
        }
    return report


def random_checks(w_values: tuple[int, ...], samples: int, seed: int) -> dict[str, object]:
    rng = random.Random(seed)
    report: dict[str, object] = {}
    for w in w_values:
        universe = tuple(range(1, 1 << w))
        diameter_histogram: Counter[int] = Counter()
        geodesic_count = 0
        accepted = 0
        attempts = 0
        while accepted < samples:
            attempts += 1
            probability = rng.uniform(0.08, 0.65)
            support = tuple(x for x in universe if rng.random() < probability)
            if rank_binary(support, w) != w:
                continue
            diameter, checked = check_support(support, w)
            accepted += 1
            geodesic_count += checked
            diameter_histogram[diameter] += 1
        report[str(w)] = {
            "samples": accepted,
            "sampling_attempts": attempts,
            "diametral_geodesics_checked": geodesic_count,
            "diameter_histogram": dict(sorted(diameter_histogram.items())),
        }
    return report


def endpoint_example() -> dict[str, object]:
    # In F_2^3 take B=(e1,e2) and the complete affine fibre e3+span(B).
    # The Cayley diameter is two and that fibre has four points.
    support = (1, 2, 4, 5, 6, 7)
    distances = word_lengths(support, 3)
    fibre = (4, 5, 6, 7)
    return {
        "w": 3,
        "support": support,
        "diameter": max(distances),
        "geodesic_basis": (1, 2),
        "nonzero_fibre": fibre,
        "fibre_size": len(fibre),
        "D_plus_one": max(distances) + 1,
        "pairwise_coordinate_distances": sorted(
            {bin(x ^ y).count("1") for x, y in itertools.combinations(fibre, 2)}
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-exhaustive-w", type=int, default=4)
    parser.add_argument("--random-samples", type=int, default=250)
    parser.add_argument("--seed", type=int, default=20260816)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name(
            "phase3_geodesic_fibre_bound_results.json"
        ),
    )
    args = parser.parse_args()

    result = {
        "claim": (
            "Every diametral geodesic has an independent generator set; its "
            "span contains no other generators; all nonzero affine fibres "
            "have coordinate diameter at most two; for D>=3 they have size "
            "at most D+1."
        ),
        "exhaustive": exhaustive(args.max_exhaustive_w),
        "random": random_checks((5, 6), args.random_samples, args.seed),
        "diameter_two_endpoint": endpoint_example(),
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
