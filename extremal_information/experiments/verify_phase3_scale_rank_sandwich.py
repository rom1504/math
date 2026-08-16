#!/usr/bin/env python3
"""Finite audit of the scale-rank sandwich and its sharp thresholds."""

from __future__ import annotations

from itertools import product
import json
from pathlib import Path


def add(x, y, q):
    return tuple((a + b) % q for a, b in zip(x, y))


def scale(a, x, q):
    return tuple((a * b) % q for b in x)


def span(rows, q):
    if not rows:
        return frozenset({(0,) * 0})
    d = len(rows[0])
    return frozenset(
        sum_vectors((scale(a, row, q) for a, row in zip(coeffs, rows)), d, q)
        for coeffs in product(range(q), repeat=len(rows))
    )


def sum_vectors(vectors, d, q):
    out = (0,) * d
    for vector in vectors:
        out = add(out, vector, q)
    return out


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


def all_subspaces_binary(d):
    words = tuple(product(range(2), repeat=d))
    zero = (0,) * d
    subspaces = {frozenset({zero})}
    frontier = [frozenset({zero})]
    while frontier:
        current = frontier.pop()
        for vector in words:
            if vector in current:
                continue
            enlarged = frozenset(current | {add(x, vector, 2) for x in current})
            if enlarged not in subspaces:
                subspaces.add(enlarged)
                frontier.append(enlarged)
    return tuple(subspaces)


def dimension(subspace, q):
    size = len(subspace)
    out = 0
    while q**out < size:
        out += 1
    assert q**out == size
    return out


def min_nonzero_norm(subspace, norm):
    values = [norm(x) for x in subspace if any(x)]
    return min(values) if values else float("inf")


def two_scale_rank_curve():
    q, d, r, level = 2, 5, 2, 7
    subspaces = all_subspaces_binary(d)

    def norm(x):
        return level * int(any(x[:r])) + int(any(x))

    thresholds = (0.5, 1, level, level + 1)
    observed = {}
    expected = {0.5: d, 1: r, level: r, level + 1: 0}
    for threshold in thresholds:
        value = max(
            dimension(c, q)
            for c in subspaces
            if min_nonzero_norm(c, norm) > threshold
        )
        assert value == expected[threshold]
        observed[str(threshold)] = value

    # The non-strict analogue fails: ker(pi) has minimum exactly one and
    # dimension D-r > r.
    kernel = frozenset(x for x in product(range(q), repeat=d) if not any(x[:r]))
    assert min_nonzero_norm(kernel, norm) == 1
    assert dimension(kernel, q) == d - r > r
    return len(subspaces), observed, dimension(kernel, q)


def profile(u, vectors, q, norm):
    return min(
        2 * sum(a != 0 for a in z)
        + norm(add(u, sum_vectors((scale(a, v, q) for a, v in zip(z, vectors)), len(u), q), q))
        for z in product(range(q), repeat=len(vectors))
    )


def two_scale_decoder():
    q, d, r, k, level = 2, 5, 2, 1, 7
    words = tuple(product(range(q), repeat=d))

    def norm(x):
        return level * int(any(x[:r])) + int(any(x))

    bases = [(v,) for v in words if any(v)]
    max_error = 0
    projected_states = set()
    for vectors in bases:
        projected = frozenset(
            scale(a, vectors[0][:r], q)
            for a in range(q)
        )
        projected_states.add(projected)
        for u in words:
            decoder = (level + 1) * int(u[:r] not in projected)
            error = profile(u, vectors, q, norm) - decoder
            assert 0 <= error <= 3  # a+b+2k = 1+0+2.
            max_error = max(max_error, error)

    # A section of pi contains all one-dimensional projected carriers.
    section_vectors = [v for v in words if any(v[:r]) and not any(v[r:])]
    section_profiles = {
        v: tuple(profile(u, (v,), q, norm) for u in words)
        for v in section_vectors
    }
    minimum_separation = min(
        max(abs(a - b) for a, b in zip(section_profiles[v], section_profiles[vp]))
        for i, v in enumerate(section_vectors)
        for vp in section_vectors[i + 1 :]
    )
    assert minimum_separation > level - 2 * k
    return len(bases), len(projected_states), max_error, len(section_vectors), minimum_separation


def hamming_singleton():
    q, d = 2, 4
    subspaces = all_subspaces_binary(d)

    def hamming(x):
        return sum(x)

    values = {}
    for h in range(d + 1):
        separated_rank = max(
            dimension(c, q)
            for c in subspaces
            if min_nonzero_norm(c, hamming) > h
        )
        assert separated_rank <= d - h
        values[str(h)] = separated_rank
    return values


def presentation_radius_sharpness():
    k = 3
    d = 3 * k
    vectors = []
    for j in range(k):
        vectors.append(tuple(int(3 * j <= i < 3 * j + 3) for i in range(d)))
    query = sum_vectors(vectors, d, 2)
    value = profile(query, tuple(vectors), 2, lambda x: sum(x))
    assert value == 2 * k
    return {"k": k, "profile_value": value, "distance_to_image": 0}


def additive_constant_sharpness():
    # X=F_2^2 with weighted l1 norm; Y=F_2 with unit norm.
    fibre_a, lift_b, presentation_p = 4, 3, 5

    def xnorm(x):
        return (1 + lift_b) * x[0] + fibre_a * x[1]

    query = (0, 0)
    carrier = (1, 1)
    quotient_decoder = 1
    response = xnorm(add(query, carrier, 2)) + presentation_p
    error = response - quotient_decoder
    assert error == fibre_a + lift_b + presentation_p
    return {
        "fibre_term": fibre_a,
        "lift_term": lift_b,
        "presentation_term": presentation_p,
        "attained_error": error,
    }


def main():
    subspace_count, curve, nonstrict_kernel_dimension = two_scale_rank_curve()
    maps, states, decoder_error, packed, separation = two_scale_decoder()
    result = {
        "two_scale_rank_curve": {
            "binary_subspaces_enumerated": subspace_count,
            "observed_s_W": curve,
            "nonstrict_kernel_dimension": nonstrict_kernel_dimension,
        },
        "two_scale_decoder": {
            "maps_checked": maps,
            "projected_subspace_states": states,
            "maximum_decoder_error": decoder_error,
            "proved_error": 3,
            "section_profiles_packed": packed,
            "minimum_response_separation": separation,
            "proved_strict_lower_threshold": 7 - 2,
        },
        "hamming_separated_ranks": hamming_singleton(),
        "presentation_radius_sharpness": presentation_radius_sharpness(),
        "additive_constant_sharpness": additive_constant_sharpness(),
        "status": "all assertions passed",
    }
    out = Path(__file__).with_name("phase3_scale_rank_sandwich_results.json")
    out.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
