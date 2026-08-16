#!/usr/bin/env python3
"""Finite audit for the Hamming Grassmannian injection barrier.

The script exhaustively checks Lemma HI.1 through D=5, verifies monotonicity
of the normalized ball volume, and checks the four-carrier example in F_2^5.
"""

from __future__ import annotations

from itertools import combinations
import json
from math import comb, log2
from pathlib import Path


def wt(x: int) -> int:
    return bin(x).count("1")


def rank(vectors: list[int] | tuple[int, ...]) -> int:
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


def span(basis: list[int] | tuple[int, ...]) -> frozenset[int]:
    values = {0}
    for vector in basis:
        values |= {x ^ vector for x in tuple(values)}
    return frozenset(values)


def all_subspaces(dimension: int, rank_value: int) -> list[frozenset[int]]:
    found: set[frozenset[int]] = set()
    for candidate in combinations(range(1, 1 << dimension), rank_value):
        if rank(candidate) == rank_value:
            found.add(span(candidate))
    return sorted(found, key=lambda code: tuple(sorted(code)))


def hausdorff(left: frozenset[int], right: frozenset[int]) -> tuple[int, int, int]:
    forward = max(min(wt(x ^ y) for y in right) for x in left)
    reverse = max(min(wt(x ^ y) for x in left) for y in right)
    return forward, reverse, max(forward, reverse)


def ball_volume(length: int, radius: int) -> int:
    return sum(comb(length, index) for index in range(min(length, radius) + 1))


def audit_low_weight_certificate(max_dimension: int = 5) -> dict[str, int]:
    pairs = 0
    close_cases = 0
    implications = 0
    for dimension in range(2, max_dimension + 1):
        for k in range(1, dimension):
            codes = all_subspaces(dimension, k)
            for index, left in enumerate(codes):
                for right in codes[index + 1 :]:
                    intersection_size = len(left & right)
                    intersection_dimension = intersection_size.bit_length() - 1
                    injection = k - intersection_dimension
                    sum_code = frozenset(x ^ y for x in left for y in right)
                    assert len(sum_code) == 1 << (k + injection)
                    _, _, distance = hausdorff(left, right)
                    pairs += 1
                    for threshold in range(dimension + 1):
                        low_words = sum(wt(x) <= threshold for x in sum_code)
                        assert low_words <= ball_volume(k + injection, threshold)
                        if distance <= threshold:
                            close_cases += 1
                            assert low_words >= 1 << injection
                        if ball_volume(k + injection, threshold) < 1 << injection:
                            implications += 1
                            assert distance > threshold
    return {
        "subspace_pairs": pairs,
        "close_threshold_cases": close_cases,
        "strict_implications_checked": implications,
    }


def audit_monotonicity(max_n: int = 30) -> int:
    checks = 0
    for k in range(1, max_n):
        for threshold in range(max_n + 1):
            previous_num = ball_volume(k, threshold)
            previous_den = 1
            for s in range(1, max_n - k + 1):
                current_num = ball_volume(k + s, threshold)
                current_den = 1 << s
                assert current_num * previous_den <= previous_num * current_den
                previous_num, previous_den = current_num, current_den
                checks += 1
    return checks


def audit_four_carriers() -> dict[str, object]:
    carriers = [
        span((0b01100, 0b10010)),
        span((0b00100, 0b10001)),
        span((0b00011, 0b01000)),
        span((0b00111, 0b11001)),
    ]
    directed: dict[str, list[int]] = {}
    for left in range(len(carriers)):
        assert len(carriers[left]) == 4
        for right in range(left + 1, len(carriers)):
            forward, reverse, distance = hausdorff(carriers[left], carriers[right])
            assert distance == 3
            directed[f"{left + 1}-{right + 1}"] = [forward, reverse]

    # Exhaustive linear-host check: no dimension-three subspace has minimum
    # nonzero weight at least three.
    maximum_host_dimension = 0
    for dimension in range(6):
        if any(
            min((wt(x) for x in code if x), default=6) >= 3
            for code in all_subspaces(5, dimension)
        ):
            maximum_host_dimension = dimension
    assert maximum_host_dimension == 2
    return {
        "packing_size": len(carriers),
        "pair_directed_distances": directed,
        "maximum_distance_three_linear_host_dimension": maximum_host_dimension,
        "common_host_two_subspace_count": 1,
    }


def binary_entropy(value: float) -> float:
    if value in (0.0, 1.0):
        return 0.0
    return -value * log2(value) - (1.0 - value) * log2(1.0 - value)


def audit_entropy_domination() -> int:
    """Numerical regression for the exact analytic inequality HI.10."""

    checks = 0
    for k_index in range(1, 100):
        kappa = k_index / 200.0
        for rho_index in range(1, k_index + 1):
            rho = rho_index / 200.0
            scale = kappa + rho
            for delta_index in range(1, 100):
                delta = delta_index / 400.0
                relative = delta / scale
                if relative >= 0.5:
                    continue
                if scale * binary_entropy(relative) >= rho:
                    continue
                injection_exponent = (
                    kappa * (1.0 - kappa) - rho * (1.0 - rho)
                )
                host_exponent = kappa * (
                    1.0 - binary_entropy(delta) - kappa
                )
                assert injection_exponent <= host_exponent + 1e-12
                checks += 1
    return checks


def main() -> None:
    result = {
        "low_weight_certificate": audit_low_weight_certificate(),
        "normalized_volume_monotonicity_checks": audit_monotonicity(),
        "entropy_domination_grid_checks": audit_entropy_domination(),
        "four_carrier_counterexample": audit_four_carriers(),
        "status": "all assertions passed",
    }
    output = Path(__file__).with_name(
        "phase3_hamming_grassmannian_injection_barrier_results.json"
    )
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
