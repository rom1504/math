#!/usr/bin/env python3
"""Verify the Phase-2 syndrome-profile composition examples.

The integer ``s`` encodes a vector of F_2^w.  A parity-check fragment is
represented only by the set of its nonzero column types.  The script checks
the min-plus convolution law, the special environments that expose every
support bit, the strict-quotient code pair, and the outer-spectrum collision.
"""

from __future__ import print_function

import argparse
import json
from collections import deque
from itertools import combinations
from pathlib import Path


INF = 10 ** 6


def syndrome_profile(w, columns):
    """Shortest subset-sum length for every syndrome."""
    columns = tuple(sorted(set(columns)))
    distance = [INF] * (1 << w)
    distance[0] = 0
    queue = deque([0])
    while queue:
        x = queue.popleft()
        for column in columns:
            y = x ^ column
            if distance[y] > distance[x] + 1:
                distance[y] = distance[x] + 1
                queue.append(y)
    return tuple(distance)


def min_convolution(left, right):
    size = len(left)
    return tuple(
        min(left[u] + right[s ^ u] for u in range(size))
        for s in range(size)
    )


def histogram(profile):
    finite = [value for value in profile if value < INF]
    result = [0] * (max(finite) + 1)
    for value in finite:
        result[value] += 1
    return tuple(result)


def kernel_weight_enumerator(columns):
    """Weight enumerator of ker(H), where ``columns`` are the columns of H."""
    n = len(columns)
    result = [0] * (n + 1)
    for mask in range(1 << n):
        syndrome = 0
        for i, column in enumerate(columns):
            if (mask >> i) & 1:
                syndrome ^= column
        if syndrome == 0:
            result[bin(mask).count("1")] += 1
    return tuple(result)


def is_spanning(w, support):
    return max(syndrome_profile(w, support)) < INF


def exhaustive_small_checks():
    counts = {}
    for w in (2, 3):
        nonzero = tuple(range(1, 1 << w))
        supports = []
        for mask in range(1 << len(nonzero)):
            support = tuple(
                value for i, value in enumerate(nonzero) if (mask >> i) & 1
            )
            if is_spanning(w, support):
                supports.append(support)

        # Every special environment E_s reports whether s is in the support.
        for support in supports:
            for s in nonzero:
                environment = tuple(value for value in nonzero if value != s)
                radius = max(syndrome_profile(w, support + environment))
                assert radius == (1 if s in support else 2)

        # Concatenating fragments is min-plus convolution, and for binary
        # Hamming correction this equals taking the union of column types.
        for left in supports:
            left_profile = syndrome_profile(w, left)
            for right in supports:
                right_profile = syndrome_profile(w, right)
                assert min_convolution(left_profile, right_profile) == \
                    syndrome_profile(w, left + right)

        counts[str(w)] = {
            "spanning_supports": len(supports),
            "ordered_compositions_checked": len(supports) ** 2,
            "special_environment_checks": len(supports) * len(nonzero),
        }
    return counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "extremal_information/experiments/phase2_code_syndrome_profiles_results.json"
        ),
    )
    args = parser.parse_args()

    # Same syndrome state, but nonisometric length-five codes.
    strict_a = (1, 1, 1, 2, 3)
    strict_b = (1, 1, 2, 2, 3)
    strict_profile_a = syndrome_profile(2, strict_a)
    strict_profile_b = syndrome_profile(2, strict_b)
    assert strict_profile_a == strict_profile_b == (0, 1, 1, 1)
    strict_enum_a = kernel_weight_enumerator(strict_a)
    strict_enum_b = kernel_weight_enumerator(strict_b)
    assert strict_enum_a == (1, 0, 3, 3, 0, 1)
    assert strict_enum_b == (1, 0, 2, 4, 1, 0)

    # Same root-averaged outer spectrum, but a fixed appended fragment
    # separates the covering radii.
    outer_a = (1, 2, 3, 4)
    outer_b = (1, 2, 4, 7)
    environment = (1, 3, 5, 6)
    profile_a = syndrome_profile(3, outer_a)
    profile_b = syndrome_profile(3, outer_b)
    assert profile_a == (0, 1, 1, 1, 1, 2, 2, 2)
    assert profile_b == (0, 1, 1, 2, 1, 2, 2, 1)
    assert histogram(profile_a) == histogram(profile_b) == (1, 4, 3)
    radius_a = max(syndrome_profile(3, outer_a + environment))
    radius_b = max(syndrome_profile(3, outer_b + environment))
    assert (radius_a, radius_b) == (2, 1)

    output = {
        "exhaustive_checks": exhaustive_small_checks(),
        "strict_quotient_pair": {
            "columns_a": strict_a,
            "columns_b": strict_b,
            "common_syndrome_profile": strict_profile_a,
            "kernel_weight_enumerator_a": strict_enum_a,
            "kernel_weight_enumerator_b": strict_enum_b,
        },
        "outer_spectrum_collision": {
            "columns_a": outer_a,
            "columns_b": outer_b,
            "environment_columns": environment,
            "profile_a": profile_a,
            "profile_b": profile_b,
            "common_profile_histogram": histogram(profile_a),
            "common_outer_polynomial_coefficients": (2, 8, 6),
            "composite_covering_radii": (radius_a, radius_b),
        },
    }
    rendered = json.dumps(output, indent=2, sort_keys=True) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(rendered)
    print(rendered, end="")


if __name__ == "__main__":
    main()
