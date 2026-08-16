#!/usr/bin/env python3
"""Exhaustive small checks for two phase-3 closed-algebra theorems.

The default run verifies:

1. the exact future-response metric on every pair of flats of the complete
   binary projective matroid through width four, against every future flat;
2. the carrier-span pointwise and radius inequalities through width three,
   for one to three carriers, every background dictionary, every future
   dictionary, and every target.

Integers encode vectors in F_2^w.  A dictionary mask uses bit ``v-1`` for
the nonzero vector ``v``.  The script is deterministic and writes only its
JSON result beside itself unless ``--output`` is supplied.
"""

from __future__ import annotations

import argparse
import json
from functools import lru_cache
from itertools import combinations, combinations_with_replacement
from pathlib import Path


def span(vectors: object) -> frozenset[int]:
    values = {0}
    for vector in vectors:
        values |= {x ^ int(vector) for x in tuple(values)}
    return frozenset(values)


def all_subspaces(w: int) -> list[frozenset[int]]:
    """Enumerate all subspaces by closing all independent sets of size <=w."""

    seen: set[frozenset[int]] = {frozenset({0})}
    nonzero = range(1, 1 << w)
    for size in range(1, w + 1):
        for candidate in combinations(nonzero, size):
            candidate_span = span(candidate)
            if len(candidate_span) == 1 << size:
                seen.add(candidate_span)
    return sorted(seen, key=lambda space: (len(space), tuple(sorted(space))))


def dimension(space: frozenset[int]) -> int:
    return len(space).bit_length() - 1


def atom_mask(space: frozenset[int]) -> int:
    mask = 0
    for vector in space:
        if vector:
            mask |= 1 << (vector - 1)
    return mask


def mask_vectors(mask: int, w: int) -> tuple[int, ...]:
    return tuple(v for v in range(1, 1 << w) if mask & (1 << (v - 1)))


def check_matroid_metric(max_w: int) -> list[dict[str, int]]:
    cases: list[dict[str, int]] = []
    for w in range(1, max_w + 1):
        flats = all_subspaces(w)
        ordered_pairs = 0
        response_evaluations = 0
        for X in flats:
            rank_x = dimension(X)
            for Y in flats:
                rank_y = dimension(Y)
                join = span(X | Y)
                formula = max(
                    dimension(join) - rank_x,
                    dimension(join) - rank_y,
                )
                actual = 0
                for T in flats:
                    response_evaluations += 1
                    gap = abs(dimension(span(X | T)) - dimension(span(Y | T)))
                    actual = max(actual, gap)
                assert actual == formula, (w, X, Y, actual, formula)
                ordered_pairs += 1
        cases.append(
            {
                "w": w,
                "flat_count": len(flats),
                "ordered_pair_count": ordered_pairs,
                "future_response_evaluations": response_evaluations,
            }
        )
    return cases


def check_multi_carrier(max_w: int, max_m: int) -> list[dict[str, int]]:
    cases: list[dict[str, int]] = []
    for w in range(1, max_w + 1):
        all_atoms = (1 << ((1 << w) - 1)) - 1
        spaces = all_subspaces(w)
        carrier_masks = [atom_mask(space) for space in spaces]

        @lru_cache(maxsize=None)
        def profile(mask: int) -> tuple[int, ...]:
            generators = mask_vectors(mask, w)
            distances = [w + 1] * (1 << w)
            distances[0] = 0
            frontier = [0]
            depth = 0
            while frontier:
                next_frontier: list[int] = []
                for x in frontier:
                    for generator in generators:
                        y = x ^ generator
                        if distances[y] > depth + 1:
                            distances[y] = depth + 1
                            next_frontier.append(y)
                frontier = next_frontier
                depth += 1
            return tuple(distances)

        for m in range(1, max_m + 1):
            # Each unique pair is (D, Dbar), after allowing every background.
            source_pairs: set[tuple[int, int]] = set()
            for carriers in combinations_with_replacement(carrier_masks, m):
                dense_union = 0
                carrier_vectors: list[int] = []
                for carrier in carriers:
                    dense_union |= carrier
                    carrier_vectors.extend(mask_vectors(carrier, w))
                span_carrier = atom_mask(span(carrier_vectors))
                for background in range(all_atoms + 1):
                    D = background | dense_union
                    Dbar = background | span_carrier
                    # Work in the finite-radius dictionary model.
                    if max(profile(D)) <= w:
                        source_pairs.add((D, Dbar))

            pointwise_checks = 0
            radius_checks = 0
            for D, Dbar in source_pairs:
                for future in range(all_atoms + 1):
                    small = profile(D | future)
                    large = profile(Dbar | future)
                    for target in range(1 << w):
                        gap = small[target] - large[target]
                        assert 0 <= gap <= m - 1, (
                            w,
                            m,
                            D,
                            Dbar,
                            future,
                            target,
                            gap,
                        )
                        pointwise_checks += 1
                    radius_gap = max(small) - max(large)
                    assert 0 <= radius_gap <= m - 1
                    radius_checks += 1

            cases.append(
                {
                    "w": w,
                    "carrier_count": m,
                    "unique_spanning_source_pairs": len(source_pairs),
                    "future_dictionary_count": all_atoms + 1,
                    "pointwise_checks": pointwise_checks,
                    "radius_checks": radius_checks,
                }
            )
    return cases


def run(matroid_max_w: int, carrier_max_w: int, carrier_max_m: int) -> dict[str, object]:
    return {
        "status": "passed",
        "matroid_response_metric": check_matroid_metric(matroid_max_w),
        "multi_carrier_arbitrary_future": check_multi_carrier(
            carrier_max_w, carrier_max_m
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matroid-max-w", type=int, default=4)
    parser.add_argument("--carrier-max-w", type=int, default=3)
    parser.add_argument("--carrier-max-m", type=int, default=3)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name(
            "phase3_closed_algebra_claims_results.json"
        ),
    )
    args = parser.parse_args()
    result = run(args.matroid_max_w, args.carrier_max_w, args.carrier_max_m)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
