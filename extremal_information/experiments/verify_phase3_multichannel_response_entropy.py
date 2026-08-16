#!/usr/bin/env python3
"""Finite checks for phase3_multichannel_response_entropy.md.

No external packages are required.  All vectors are integers with xor as
addition and ``bit_count`` as Hamming weight.
"""

from __future__ import annotations

import itertools
import json
from pathlib import Path


def wt(x: int) -> int:
    return bin(x).count("1")


def profile(d: int, columns: tuple[int, ...]) -> tuple[int, ...]:
    centers: list[tuple[int, int]] = []
    for mask in range(1 << len(columns)):
        value = 0
        for j, column in enumerate(columns):
            if (mask >> j) & 1:
                value ^= column
        centers.append((2 * wt(mask), value))
    return tuple(
        min(cost + wt(u ^ value) for cost, value in centers)
        for u in range(1 << d)
    )


def span(basis: tuple[int, ...]) -> frozenset[int]:
    values = {0}
    for vector in basis:
        values |= {x ^ vector for x in tuple(values)}
    return frozenset(values)


def subspaces_of_dimension(code: frozenset[int], k: int) -> list[frozenset[int]]:
    out: set[frozenset[int]] = set()
    for basis in itertools.combinations(sorted(code - {0}), k):
        space = span(basis)
        if len(space) == 1 << k:
            out.add(space)
    return sorted(out, key=lambda c: tuple(sorted(c)))


def lex_basis(code: frozenset[int], k: int) -> tuple[int, ...]:
    chosen: list[int] = []
    current = frozenset({0})
    for value in sorted(code - {0}):
        if value not in current:
            chosen.append(value)
            current = span(tuple(chosen))
            if len(chosen) == k:
                return tuple(chosen)
    raise AssertionError("code has insufficient dimension")


def check_exact_collision_classification() -> dict[str, int]:
    stats: dict[str, int] = {}
    for d, k in ((4, 2), (5, 2)):
        seen: set[tuple[int, ...]] = set()
        expected_sets: set[frozenset[int]] = set()
        for columns in itertools.product(range(1 << d), repeat=k):
            f = profile(d, columns)
            heavy = frozenset(v for v in columns if wt(v) >= 3)
            recovered = frozenset(
                u
                for u, value in enumerate(f)
                if wt(u) >= 3 and value == 2
            )
            assert recovered == heavy
            seen.add(f)
            expected_sets.add(heavy)
        assert len(seen) == len(expected_sets)
        heavy_count = sum(1 for u in range(1 << d) if wt(u) >= 3)
        expected_count = sum(
            __import__("math").comb(heavy_count, j) for j in range(k + 1)
        )
        assert len(seen) == expected_count
        stats[f"D{d}_k{k}_profiles"] = len(seen)
    return stats


def check_basis_dependence() -> dict[str, int]:
    t = 3
    # Blocks X,Y,Z.  a is one on X,Z; b is one on Y,Z.
    x = (1 << t) - 1
    a = x | (x << (2 * t))
    b = (x << t) | (x << (2 * t))
    f_ab = profile(3 * t, (a, b))
    f_aab = profile(3 * t, (a, a ^ b))
    assert span((a, b)) == span((a, a ^ b))
    assert f_ab[b] == 2
    assert f_aab[b] == 4
    return {"D": 3 * t, "same_span_profile_gap_at_b": f_aab[b] - f_ab[b]}


def check_span_sandwich() -> dict[str, int]:
    checked = 0
    worst_upper_gap = 0
    for d in range(1, 6):
        for k in range(0, 3):
            for columns in itertools.product(range(1 << d), repeat=k):
                f = profile(d, columns)
                code = span(columns)
                rank = (len(code)).bit_length() - 1
                for u, value in enumerate(f):
                    dc = min(wt(u ^ c) for c in code)
                    assert dc <= value <= dc + 2 * rank
                    worst_upper_gap = max(worst_upper_gap, value - dc)
                    checked += 1
    return {"pointwise_checks": checked, "worst_upper_gap": worst_upper_gap}


def simplex_code(s: int, repeats: int) -> frozenset[int]:
    """Binary simplex [repeats*(2^s-1), s, repeats*2^(s-1)] code."""
    words: set[int] = set()
    nonzero_points = list(range(1, 1 << s))
    for linear_form in range(1 << s):
        word = 0
        coordinate = 0
        for _ in range(repeats):
            for point in nonzero_points:
                if wt(linear_form & point) & 1:
                    word |= 1 << coordinate
                coordinate += 1
        words.add(word)
    return frozenset(words)


def check_good_code_subspace_packing() -> dict[str, object]:
    s, repeats = 3, 2
    d = repeats * ((1 << s) - 1)
    code0 = simplex_code(s, repeats)
    minimum_distance = min(wt(x) for x in code0 - {0})
    results: dict[str, object] = {
        "D": d,
        "host_dimension": s,
        "host_minimum_distance": minimum_distance,
    }
    for k in (1, 2):
        codes = subspaces_of_dimension(code0, k)
        profiles = [(code, profile(d, lex_basis(code, k))) for code in codes]
        min_gap = d
        min_hausdorff = d
        for i, (code, f) in enumerate(profiles):
            for other, g in profiles[i + 1 :]:
                gap = max(abs(a - b) for a, b in zip(f, g))
                hausdorff = max(
                    max(min(wt(x ^ y) for y in other) for x in code),
                    max(min(wt(x ^ y) for x in code) for y in other),
                )
                min_gap = min(min_gap, gap)
                min_hausdorff = min(min_hausdorff, hausdorff)
                assert gap >= minimum_distance - 2 * k
                assert hausdorff >= minimum_distance
        results[f"k{k}"] = {
            "subspaces": len(codes),
            "minimum_profile_gap": min_gap,
            "minimum_hausdorff_distance": min_hausdorff,
            "proved_gap_floor": minimum_distance - 2 * k,
        }
    return results


def check_two_fragment_identity() -> dict[str, int]:
    # Directly enumerate selections from paired quotient columns.  A mask on
    # 2k letters has quotient zero iff its two bits agree in each channel.
    checks = 0
    for d in range(1, 5):
        for k in range(1, 4):
            for columns in itertools.product(range(1 << d), repeat=k):
                f = profile(d, columns)
                for u in range(1 << d):
                    best = d
                    for mask in range(1 << (2 * k)):
                        quotient = 0
                        kernel = 0
                        used = wt(mask)
                        for j, v in enumerate(columns):
                            if (mask >> (2 * j)) & 1:
                                quotient ^= 1 << j
                            if (mask >> (2 * j + 1)) & 1:
                                quotient ^= 1 << j
                                kernel ^= v
                        if quotient == 0:
                            best = min(best, used + wt(u ^ kernel))
                    assert best == f[u]
                    checks += 1
    return {"endpoint_checks": checks}


def main() -> None:
    result = {
        "exact_collision_classification": check_exact_collision_classification(),
        "basis_dependence": check_basis_dependence(),
        "span_sandwich": check_span_sandwich(),
        "good_code_subspace_packing": check_good_code_subspace_packing(),
        "two_fragment_identity": check_two_fragment_identity(),
    }
    output = Path(__file__).with_name("phase3_multichannel_response_entropy_results.json")
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
