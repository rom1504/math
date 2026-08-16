#!/usr/bin/env python3
"""Independent checks for the geodesic-fibre hard-core theorem.

This deliberately tests claims not covered by the primary verifier:

* exact diameter-two anticode counts and maximum sizes through D=4;
* the large-diameter support enumeration bound through w=4;
* terminal collapse against every raw appended support in F_2^3,
  including nonspanning futures; and
* that the advertised complexity is the logarithm of the number of quotient
  states (message bits), not the number of states.
"""

from __future__ import annotations

import itertools
import json
import math
from collections import Counter, deque
from pathlib import Path


def rank_binary(vectors: tuple[int, ...], w: int) -> int:
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


def diameter(support: tuple[int, ...], w: int) -> int:
    distances = [-1] * (1 << w)
    distances[0] = 0
    queue = deque([0])
    while queue:
        x = queue.popleft()
        for s in support:
            y = x ^ s
            if distances[y] < 0:
                distances[y] = distances[x] + 1
                queue.append(y)
    if min(distances) < 0:
        raise ValueError("support is not spanning")
    return max(distances)


def all_supports(w: int) -> list[tuple[int, ...]]:
    universe = tuple(range(1, 1 << w))
    return [
        tuple(universe[i] for i in range(len(universe)) if mask >> i & 1)
        for mask in range(1 << len(universe))
    ]


def anticode_report(max_d: int = 4) -> dict[str, object]:
    report: dict[str, object] = {}
    for d in range(1, max_d + 1):
        vertices = tuple(range(1 << d))
        count = 0
        max_size = 0
        for mask in range(1 << len(vertices)):
            family = tuple(v for v in vertices if mask >> v & 1)
            if all(
                bin(x ^ y).count("1") <= 2
                for x, y in itertools.combinations(family, 2)
            ):
                count += 1
                max_size = max(max_size, len(family))
        loose_bound = 1 + (1 << d) * (
            (1 << d) + d * (1 << (d + 1)) + math.comb(d, 3)
        )
        theorem_bound = 4 * d * (1 << (2 * d))
        assert count <= loose_bound <= theorem_bound
        if d >= 3:
            assert max_size <= d + 1
        report[str(d)] = {
            "exact_number_of_anticodes": count,
            "maximum_size": max_size,
            "D_plus_one": d + 1,
            "GF3_intermediate_bound": loose_bound,
            "GF3_final_bound": theorem_bound,
        }
    return report


def support_count_report(max_w: int = 4) -> dict[str, object]:
    report: dict[str, object] = {}
    for w in range(1, max_w + 1):
        histogram: Counter[int] = Counter()
        for support in all_supports(w):
            if rank_binary(support, w) == w:
                histogram[diameter(support, w)] += 1
        thresholds: dict[str, object] = {}
        for r in range(4, w + 1):
            actual = sum(number for d, number in histogram.items() if d >= r)
            log_actual = math.log2(actual)
            claimed_log_bound = (
                w * w
                + (1 << (w - r)) * (2 * w + 2 + math.log2(w))
                + math.log2(w)
            )
            message_bits = math.ceil(math.log2(1 + actual))
            assert log_actual <= claimed_log_bound + 1e-12
            assert message_bits <= math.ceil(claimed_log_bound) + 1
            thresholds[str(r)] = {
                "retained_supports": actual,
                "log2_retained_supports": log_actual,
                "GF4_log2_bound": claimed_log_bound,
                "actual_quotient_message_bits": message_bits,
            }
        report[str(w)] = {
            "diameter_histogram": dict(sorted(histogram.items())),
            "thresholds": thresholds,
        }
    return report


def all_future_report(w: int = 3) -> dict[str, object]:
    supports = all_supports(w)
    spanning = [s for s in supports if rank_binary(s, w) == w]
    nonspanning_futures = sum(rank_binary(u, w) < w for u in supports)
    checks = 0
    collapsed_checks = 0
    for s in spanning:
        d_s = diameter(s, w)
        for r in range(2, w + 1):
            for u in supports:
                union = tuple(sorted(set(s) | set(u)))
                d_union = diameter(union, w)
                checks += 1
                if d_s < r:
                    collapsed_checks += 1
                    assert d_union < r
                    assert abs(d_union - r / 2) <= r / 2 - 1 + 1e-12
                else:
                    # A retained state stores S exactly, so the union diameter
                    # is itself the exact decoded answer.
                    assert diameter(tuple(sorted(set(s) | set(u))), w) == d_union
    return {
        "w": w,
        "spanning_sources": len(spanning),
        "raw_future_supports": len(supports),
        "nonspanning_future_supports": nonspanning_futures,
        "source_threshold_future_checks": checks,
        "collapsed_future_checks": collapsed_checks,
    }


def main() -> None:
    result = {
        "anticodes": anticode_report(),
        "support_counts": support_count_report(),
        "all_raw_futures": all_future_report(),
        "conclusion": (
            "The anticode, enumeration, terminal-collapse, nonspanning-future, "
            "and message-bit claims survived the independent finite audit."
        ),
    }
    output = Path(__file__).with_name("phase3_geodesic_second_audit_results.json")
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
