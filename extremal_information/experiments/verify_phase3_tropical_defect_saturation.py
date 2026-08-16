#!/usr/bin/env python3
"""Exact finite falsification checks for tropical defect saturation.

The checks cover three logically separate claims from the phase-three draft:

1. finite-valued min-plus kernels on small cyclic groups, including the sharp
   one-blur versus many-blur defect and finite Kleene stabilization;
2. indicator/Hamming-ball powers, with both the raw smoothing error and the
   additional one-blur-to-many-blur error recorded separately; and
3. the fixed-chart syndrome quotient for every support pair through width
   three and a deterministic sample at width four.

The program also records the excluded identically-infinite profile: convolving
it with a finite kernel leaves it identically infinite, so the theorem's norm
statement must be restricted to proper extended profiles.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
import random
from collections import deque
from pathlib import Path


INF = math.inf


def minplus_cyclic(f: tuple[float, ...], g: tuple[float, ...]) -> tuple[float, ...]:
    q = len(f)
    return tuple(min(f[y] + g[(x - y) % q] for y in range(q)) for x in range(q))


def sup_distance(f: tuple[float, ...], g: tuple[float, ...]) -> float:
    values: list[float] = []
    for left, right in zip(f, g):
        if math.isinf(left) and math.isinf(right):
            values.append(0.0)
        elif math.isinf(left) or math.isinf(right):
            return INF
        else:
            values.append(abs(left - right))
    return max(values, default=0.0)


def kernel_checks(max_q: int, max_cost: int) -> list[dict[str, int]]:
    results: list[dict[str, int]] = []
    for q in range(2, max_q + 1):
        kernel_count = 0
        profile_checks = 0
        for tail in itertools.product(range(max_cost + 1), repeat=q - 1):
            b = (0.0,) + tuple(float(x) for x in tail)
            powers = [b]
            for _ in range(1, 2 * q + 2):
                powers.append(minplus_cyclic(powers[-1], b))

            # Zero-cost padding makes powers decrease.  A shortest path on q
            # vertices is simple, so powers have stabilized by q steps.
            assert all(
                all(right <= left for left, right in zip(powers[m], powers[m + 1]))
                for m in range(len(powers) - 1)
            )
            assert all(powers[m] == powers[q] for m in range(q, len(powers)))
            b_star = powers[q]
            assert minplus_cyclic(b_star, b_star) == b_star
            assert all(b_star[x] <= b[x] for x in range(q))
            assert all(
                b_star[x] <= b_star[y] + b_star[(x - y) % q]
                for x in range(q)
                for y in range(q)
            )

            delta = sup_distance(b, b_star)
            delta_zero = (0.0,) + (INF,) * (q - 1)
            for m, b_power in enumerate(powers, start=1):
                exact_defect = sup_distance(b_power, b)
                # The delta profile attains the operator norm exactly.
                assert sup_distance(
                    minplus_cyclic(delta_zero, b_power),
                    minplus_cyclic(delta_zero, b),
                ) == exact_defect

                # Exhaust all proper profiles with entries in {0,1,2,infinity}.
                for values in itertools.product((0.0, 1.0, 2.0, INF), repeat=q):
                    if all(math.isinf(value) for value in values):
                        continue
                    left = minplus_cyclic(values, b_power)
                    right = minplus_cyclic(values, b)
                    assert sup_distance(left, right) <= exact_defect
                    profile_checks += 1
            assert max(sup_distance(power, b) for power in powers) == delta
            kernel_count += 1

        results.append(
            {
                "group_order": q,
                "kernel_count": kernel_count,
                "proper_profile_inequality_checks": profile_checks,
            }
        )
    return results


def hamming_ball_checks(max_w: int) -> list[dict[str, int]]:
    results: list[dict[str, int]] = []
    for w in range(1, max_w + 1):
        checks = 0
        for r in range(w + 1):
            for m in range(1, w + 3):
                s = min(m * r, w)
                raw_error = max(
                    weight - max(weight - s, 0) for weight in range(w + 1)
                )
                one_to_many_error = max(
                    abs(max(weight - r, 0) - max(weight - s, 0))
                    for weight in range(w + 1)
                )
                assert raw_error == s
                assert one_to_many_error == s - r
                checks += 1
        results.append({"width": w, "radius_power_checks": checks})
    return results


def word_profile(w: int, support: frozenset[int]) -> tuple[int, ...]:
    distances = [w + 1] * (1 << w)
    distances[0] = 0
    queue: deque[int] = deque([0])
    while queue:
        x = queue.popleft()
        for atom in support:
            y = x ^ atom
            if distances[y] > distances[x] + 1:
                distances[y] = distances[x] + 1
                queue.append(y)
    assert max(distances) <= w
    return tuple(distances)


def xor_convolution(f: tuple[int, ...], g: tuple[int, ...]) -> tuple[int, ...]:
    size = len(f)
    return tuple(min(f[y] + g[x ^ y] for y in range(size)) for x in range(size))


def quotient_profile(profile: tuple[int, ...], w: int, r: int) -> tuple[int, ...]:
    return tuple(
        min(profile[(coset << r) ^ h] for h in range(1 << r))
        for coset in range(1 << (w - r))
    )


def supports_with_basis(w: int) -> list[frozenset[int]]:
    basis = frozenset(1 << j for j in range(w))
    optional = [x for x in range(1, 1 << w) if x not in basis]
    return [
        basis | frozenset(optional[j] for j in range(len(optional)) if mask >> j & 1)
        for mask in range(1 << len(optional))
    ]


def fixed_chart_checks(
    exhaustive_max_w: int, sampled_w: int, sampled_pairs: int, seed: int
) -> list[dict[str, int]]:
    rng = random.Random(seed)
    results: list[dict[str, int]] = []
    for w in range(1, sampled_w + 1):
        supports = supports_with_basis(w)
        profiles = {support: word_profile(w, support) for support in supports}
        if w <= exhaustive_max_w:
            pairs = list(itertools.product(supports, repeat=2))
        else:
            pairs = [
                (rng.choice(supports), rng.choice(supports))
                for _ in range(sampled_pairs)
            ]

        quotient_checks = 0
        radius_checks = 0
        strict_collision_found = False
        for r in range(w + 1):
            cells: dict[tuple[int, ...], frozenset[int]] = {}
            for support in supports:
                coarse = quotient_profile(profiles[support], w, r)
                previous = cells.get(coarse)
                if previous is not None and previous != support:
                    strict_collision_found = True
                cells[coarse] = support

            for left, right in pairs:
                union = left | right
                union_profile = word_profile(w, union)
                coarse_union = quotient_profile(union_profile, w, r)
                coarse_product = xor_convolution(
                    quotient_profile(profiles[left], w, r),
                    quotient_profile(profiles[right], w, r),
                )
                assert coarse_union == coarse_product
                lower = max(coarse_union)
                radius = max(union_profile)
                assert lower <= radius <= lower + r
                quotient_checks += 1
                radius_checks += 1

        results.append(
            {
                "width": w,
                "support_count": len(supports),
                "ordered_support_pairs_per_chart": len(pairs),
                "chart_count": w + 1,
                "quotient_identity_checks": quotient_checks,
                "radius_interval_checks": radius_checks,
                "strict_collision_found": int(strict_collision_found),
            }
        )
    return results


def run(args: argparse.Namespace) -> dict[str, object]:
    all_infinite = (INF, INF, INF)
    finite_kernel = (0.0, 1.0, 2.0)
    excluded_output = minplus_cyclic(all_infinite, finite_kernel)
    assert all(math.isinf(value) for value in excluded_output)
    return {
        "status": "passed",
        "proper_profile_scope": {
            "identically_infinite_remains_identically_infinite": True
        },
        "finite_kernel_defect": kernel_checks(args.max_q, args.max_cost),
        "hamming_indicator": hamming_ball_checks(args.max_hamming_w),
        "fixed_chart_syndrome": fixed_chart_checks(
            args.exhaustive_chart_w,
            args.sampled_chart_w,
            args.sampled_pairs,
            args.seed,
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-q", type=int, default=4)
    parser.add_argument("--max-cost", type=int, default=3)
    parser.add_argument("--max-hamming-w", type=int, default=8)
    parser.add_argument("--exhaustive-chart-w", type=int, default=3)
    parser.add_argument("--sampled-chart-w", type=int, default=4)
    parser.add_argument("--sampled-pairs", type=int, default=512)
    parser.add_argument("--seed", type=int, default=20260816)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name(
            "phase3_tropical_defect_saturation_results.json"
        ),
    )
    args = parser.parse_args()
    result = run(args)
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
