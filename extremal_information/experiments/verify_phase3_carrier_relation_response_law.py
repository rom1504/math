#!/usr/bin/env python3
"""Finite audits for the carrier/relation response law."""

from __future__ import annotations

import heapq
from itertools import product
import json
from pathlib import Path


def add(x, y, q):
    return tuple((a + b) % q for a, b in zip(x, y))


def sub(x, y, q):
    return tuple((a - b) % q for a, b in zip(x, y))


def scale(a, x, q):
    return tuple((a * b) % q for b in x)


def rank(rows, q):
    rows = [list(row) for row in rows]
    if not rows:
        return 0
    out = 0
    for col in range(len(rows[0])):
        pivot = next((i for i in range(out, len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        rows[out], rows[pivot] = rows[pivot], rows[out]
        inv = pow(rows[out][col], -1, q)
        rows[out] = [(inv * x) % q for x in rows[out]]
        for i in range(len(rows)):
            if i != out and rows[i][col]:
                a = rows[i][col]
                rows[i] = [(x - a * y) % q for x, y in zip(rows[i], rows[out])]
        out += 1
    return out


def hamming_weight(x):
    return sum(a != 0 for a in x)


def matvec(vectors, z, q):
    out = (0,) * len(vectors[0])
    for a, v in zip(z, vectors):
        out = add(out, scale(a, v, q), q)
    return out


def code(vectors, q):
    return {
        matvec(vectors, z, q)
        for z in product(range(q), repeat=len(vectors))
    }


def profile(u, vectors, q, norm):
    return min(
        2 * hamming_weight(z) + norm(add(u, matvec(vectors, z, q), q))
        for z in product(range(q), repeat=len(vectors))
    )


def weighted_union_distances(vectors, q, norm):
    """Full rooted metric for P union R_V, for a complete weighted kernel."""
    d = len(vectors[0])
    k = len(vectors)
    zero_w = (0,) * d
    zero_q = (0,) * k
    generators = []
    for w in product(range(q), repeat=d):
        if any(w):
            generators.append((w + zero_q, norm(w)))
    for j in range(k):
        ej = tuple(int(i == j) for i in range(k))
        for a in range(1, q):
            aq = scale(a, ej, q)
            generators.append((zero_w + aq, 1))
            generators.append((scale(a, vectors[j], q) + aq, 1))

    origin = (0,) * (d + k)
    distances = {origin: 0}
    heap = [(0, origin)]
    while heap:
        current, x = heapq.heappop(heap)
        if current != distances[x]:
            continue
        for g, cost in generators:
            y = add(x, g, q)
            candidate = current + cost
            if candidate < distances.get(y, float("inf")):
                distances[y] = candidate
                heapq.heappush(heap, (candidate, y))
    return distances


def dist(u, v, q, norm):
    return norm(sub(u, v, q))


def distance_to_set(u, c, q, norm):
    return min(dist(u, v, q, norm) for v in c)


def hausdorff(c, cp, q, norm):
    return max(
        max(distance_to_set(x, cp, q, norm) for x in c),
        max(distance_to_set(x, c, q, norm) for x in cp),
    )


def check_rough_geometry(words, bases, q, norm):
    profiles = {
        vectors: tuple(profile(u, vectors, q, norm) for u in words)
        for vectors in bases
    }
    codes = {vectors: code(vectors, q) for vectors in bases}
    maximum_error = 0
    pairs = 0
    for i, vectors in enumerate(bases):
        c = codes[vectors]
        for vp in bases[i + 1 :]:
            cp = codes[vp]
            response_distance = max(
                abs(a - b) for a, b in zip(profiles[vectors], profiles[vp])
            )
            error = abs(response_distance - hausdorff(c, cp, q, norm))
            assert error <= 2 * len(vectors)
            maximum_error = max(maximum_error, error)
            pairs += 1
    return pairs, maximum_error


def discrete_collapse():
    q, d, k = 3, 3, 2
    words = list(product(range(q), repeat=d))
    bases = [v for v in product(words, repeat=k) if rank(v, q) == k]
    norm = lambda x: int(any(x))
    checks = 0
    for vectors in bases:
        for u in words:
            assert profile(u, vectors, q, norm) == int(any(u))
            checks += 1
    return len(bases), checks


def two_scale_collapse():
    q, d, k, coarse_rank, scale_size = 3, 3, 1, 1, 20
    words = list(product(range(q), repeat=d))
    bases = [(v,) for v in words if any(v)]

    def projection(x):
        return x[:coarse_rank]

    def norm(x):
        return scale_size * int(any(projection(x))) + int(any(x))

    coarse_states = set()
    maximum_error = 0
    for vectors in bases:
        c = code(vectors, q)
        projected_code = {projection(x) for x in c}
        coarse_states.add(frozenset(projected_code))
        for u in words:
            canonical = (scale_size + 1) * int(projection(u) not in projected_code)
            error = abs(profile(u, vectors, q, norm) - canonical)
            assert error <= 2 * k + 1
            maximum_error = max(maximum_error, error)
    assert len(coarse_states) == 2
    # Independently check the all-root comparison for maps with the same
    # coarse linear projection, in a smaller complete weighted Cayley graph.
    full_d = 2
    full_words = list(product(range(q), repeat=full_d))

    def full_norm(x):
        return scale_size * int(x[0] != 0) + int(any(x))

    full_metrics = {
        v: weighted_union_distances((v,), q, full_norm)
        for v in full_words
    }
    same_projection_pairs = 0
    max_full_root_gap = 0
    for i, v in enumerate(full_words):
        for vp in full_words[i + 1 :]:
            if v[0] != vp[0]:
                continue
            gap = max(abs(full_metrics[v][x] - full_metrics[vp][x]) for x in full_metrics[v])
            assert gap <= 1
            max_full_root_gap = max(max_full_root_gap, gap)
            same_projection_pairs += 1

    return (
        len(bases),
        len(coarse_states),
        maximum_error,
        same_projection_pairs,
        max_full_root_gap,
    )


def lee_geometry():
    q, d, k = 5, 2, 1
    words = list(product(range(q), repeat=d))
    # One representative for each projective line.
    seen = set()
    bases = []
    for v in words:
        if not any(v):
            continue
        c = frozenset(scale(a, v, q) for a in range(q))
        if c not in seen:
            seen.add(c)
            bases.append((v,))

    def lee(x):
        return sum(min(a, q - a) for a in x)

    pairs, maximum_error = check_rough_geometry(words, bases, q, lee)
    return len(bases), pairs, maximum_error


def flag_host():
    q, d, r, k = 2, 6, 3, 1
    words = list(product(range(q), repeat=d))

    def flag(x):
        return max((i + 1 for i, a in enumerate(x) if a), default=0)

    host_words = [x for x in words if not any(x[: d - r])]
    host_nonzero = [x for x in host_words if any(x)]
    minimum_host_weight = min(flag(x) for x in host_nonzero)
    assert minimum_host_weight == d - r + 1
    bases = [(v,) for v in host_nonzero]  # Binary lines have one basis each.
    profiles = {
        v: tuple(profile(u, v, q, flag) for u in words)
        for v in bases
    }
    minimum_response_separation = min(
        max(abs(a - b) for a, b in zip(profiles[v], profiles[vp]))
        for i, v in enumerate(bases)
        for vp in bases[i + 1 :]
    )
    assert minimum_response_separation >= minimum_host_weight - 2 * k
    pairs, maximum_error = check_rough_geometry(words, bases, q, flag)
    return (
        len(bases),
        pairs,
        minimum_host_weight,
        minimum_response_separation,
        maximum_error,
    )


def main():
    discrete_bases, discrete_checks = discrete_collapse()
    (
        two_scale_bases,
        coarse_states,
        two_scale_error,
        same_projection_pairs,
        max_full_root_gap,
    ) = two_scale_collapse()
    lee_lines, lee_pairs, lee_error = lee_geometry()
    flag_lines, flag_pairs, flag_distance, flag_separation, flag_error = flag_host()
    result = {
        "discrete_metric": {
            "field": 3,
            "D": 3,
            "k": 2,
            "independent_tuples": discrete_bases,
            "profile_values_checked": discrete_checks,
            "distinct_profiles": 1,
        },
        "two_scale_metric": {
            "field": 3,
            "D": 3,
            "k": 1,
            "maps_checked": two_scale_bases,
            "coarse_decoder_states": coarse_states,
            "max_decoder_error": two_scale_error,
            "proved_bound": 2 * 1 + 1,
            "same_projection_map_pairs_checked": same_projection_pairs,
            "max_full_root_metric_gap": max_full_root_gap,
            "proved_full_root_bound": 1,
        },
        "lee_metric": {
            "field": 5,
            "D": 2,
            "projective_lines": lee_lines,
            "line_pairs": lee_pairs,
            "max_rough_isometry_error": lee_error,
            "proved_bound": 2,
        },
        "flag_ultrametric_host": {
            "field": 2,
            "D": 6,
            "host_dimension": 3,
            "projective_lines": flag_lines,
            "line_pairs": flag_pairs,
            "host_minimum_distance": flag_distance,
            "minimum_response_separation": flag_separation,
            "max_rough_isometry_error": flag_error,
            "proved_rough_bound": 2,
        },
        "status": "all assertions passed",
    }
    out = Path(__file__).with_name("phase3_carrier_relation_response_law_results.json")
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
