#!/usr/bin/env python3
"""Finite falsifiers for benchmark_maxcut_boundary_response.md.

The script uses only exhaustive enumeration and deterministic pseudorandom
samples.  It prints a JSON report and writes no files.
"""

from __future__ import annotations

from itertools import combinations, product
import json
import math
import random


def bits(mask: int, n: int) -> tuple[int, ...]:
    return tuple((mask >> i) & 1 for i in range(n))


def hamming(x: tuple[int, ...], y: tuple[int, ...]) -> int:
    return sum(a != b for a, b in zip(x, y))


def osc(values: tuple[float, ...] | list[float]) -> float:
    return max(values) - min(values)


def response_maxcut(
    w: int, n_private: int, edges: tuple[tuple[int, int, int], ...]
) -> tuple[int, ...]:
    out = []
    for x in product((0, 1), repeat=w):
        best = -1
        for y in product((0, 1), repeat=n_private):
            z = x + y
            value = sum(weight for u, v, weight in edges if z[u] != z[v])
            best = max(best, value)
        out.append(best)
    return tuple(out)


def graph_from_mask(n: int, mask: int) -> tuple[tuple[int, int, int], ...]:
    possible = tuple(combinations(range(n), 2))
    return tuple((u, v, 1) for i, (u, v) in enumerate(possible) if mask >> i & 1)


def mapped_private_edges(
    edges: tuple[tuple[int, int, int], ...], w: int, shift: int
) -> tuple[tuple[int, int, int], ...]:
    def move(v: int) -> int:
        return v if v < w else v + shift

    return tuple((move(u), move(v), weight) for u, v, weight in edges)


def maxcut_checks() -> dict[str, object]:
    w = 3
    n_private = 2
    n = w + n_private
    edge_count = math.comb(n, 2)
    graphs = [graph_from_mask(n, mask) for mask in range(1 << edge_count)]
    responses = [response_maxcut(w, n_private, graph) for graph in graphs]

    symmetry_ok = True
    for response in responses:
        for mask in range(1 << w):
            complement = ((1 << w) - 1) ^ mask
            symmetry_ok &= response[mask] == response[complement]

    rng = random.Random(20260816)
    gluing_trials = 256
    gluing_ok = True
    for _ in range(gluing_trials):
        i = rng.randrange(len(graphs))
        j = rng.randrange(len(graphs))
        left = graphs[i]
        right = mapped_private_edges(graphs[j], w, n_private)
        direct = response_maxcut(w, 2 * n_private, left + right)
        predicted = tuple(a + b for a, b in zip(responses[i], responses[j]))
        gluing_ok &= direct == predicted

    exposure_ok = True
    exposure_tests = 0
    all_mask = (1 << w) - 1
    orbit_representatives = [mask for mask in range(1 << w) if mask <= (all_mask ^ mask)]
    for response in responses:
        M = int(osc(response)) + 1
        for target in orbit_representatives:
            target_bits = bits(target, w)
            continuation = []
            for mask in range(1 << w):
                x = bits(mask, w)
                d = min(hamming(x, target_bits), hamming(x, tuple(1 - a for a in target_bits)))
                continuation.append(M * (w - d))
            values = [a + b for a, b in zip(response, continuation)]
            maximizers = {i for i, value in enumerate(values) if value == max(values)}
            desired = {target, all_mask ^ target}
            exposure_ok &= maximizers == desired
            exposure_tests += 1

    return {
        "boundary_width": w,
        "private_vertices": n_private,
        "graphs_enumerated": len(graphs),
        "distinct_responses": len(set(responses)),
        "global_flip_symmetry": symmetry_ok,
        "gluing_trials": gluing_trials,
        "pointwise_gluing": gluing_ok,
        "orbit_exposure_tests": exposure_tests,
        "orbit_exposure": exposure_ok,
    }


def selector_response(f: tuple[int, ...], w: int) -> tuple[int, ...]:
    state_count = 1 << w
    W = max(f) - min(f)
    shifted = tuple(value - min(f) for value in f)
    L = W + 1
    C = L * w + 1
    P = 2 * (C + W) + 1
    out = []
    assignments = [bits(mask, w) for mask in range(state_count)]
    for x in assignments:
        best = -10**18
        for selector_mask in range(1 << state_count):
            selected = [a for a in range(state_count) if selector_mask >> a & 1]
            score = -C
            for a in selected:
                score += C + shifted[a] - L * hamming(x, assignments[a])
            score -= P * len(selected) * (len(selected) - 1) // 2
            best = max(best, score)
        out.append(best + min(f))
    return tuple(out)


def pairwise_and_metric_checks() -> dict[str, object]:
    rng = random.Random(1729)
    w = 3
    state_count = 1 << w
    tables = [tuple(rng.randrange(0, w + 1) for _ in range(state_count)) for _ in range(64)]
    selector_ok = all(selector_response(table, w) == table for table in tables)

    metric_ok = True
    shape_ok = True
    for f, g in zip(tables[::2], tables[1::2]):
        differences = [a - b for a, b in zip(f, g)]
        sup_distance = max(abs(value) for value in differences)
        target = max(range(state_count), key=lambda i: abs(differences[i]))
        M = max(osc(f), osc(g)) + 1
        pin = [-M * hamming(bits(i, w), bits(target, w)) for i in range(state_count)]
        exposed = abs(max(a + q for a, q in zip(f, pin)) - max(b + q for b, q in zip(g, pin)))
        metric_ok &= exposed == sup_distance

        midpoint = (max(differences) + min(differences)) / 2
        calibrated = max(abs(value - midpoint) for value in differences)
        shape_ok &= calibrated == osc(differences) / 2

    # The anchored {0,1}^(N-1) family is a 1/2-packing in shape distance.
    packing = [(0,) + tuple(v) for v in product((0, 1), repeat=state_count - 1)]
    minimum_shape_distance = math.inf
    for i, f in enumerate(packing):
        for g in packing[i + 1 :]:
            difference = [a - b for a, b in zip(f, g)]
            minimum_shape_distance = min(minimum_shape_distance, osc(difference) / 2)

    return {
        "selector_tables_checked": len(tables),
        "selector_realization": selector_ok,
        "context_metric_trials": len(tables) // 2,
        "context_sup_metric": metric_ok,
        "shape_metric": shape_ok,
        "anchored_packing_size": len(packing),
        "anchored_packing_min_distance": minimum_shape_distance,
    }


def greedy_hamming_cover(w: int, radius: int) -> list[tuple[int, ...]]:
    points = [bits(mask, w) for mask in range(1 << w)]
    uncovered = set(points)
    centres = []
    while uncovered:
        centre = max(points, key=lambda x: sum(hamming(x, y) <= radius for y in uncovered))
        centres.append(centre)
        uncovered = {y for y in uncovered if hamming(centre, y) > radius}
    return centres


def lipschitz_check() -> dict[str, object]:
    w = 6
    delta = 4.0
    radius = int(delta // 4)
    eta = delta / 2
    points = [bits(mask, w) for mask in range(1 << w)]
    anchor = points[0]
    source = points[-1]
    f = {x: float(hamming(x, source) - hamming(anchor, source)) for x in points}
    centres = greedy_hamming_cover(w, radius)
    if anchor not in centres:
        centres.append(anchor)
    rounded = {x: eta * math.ceil(f[x] / eta) for x in centres}
    g = {x: min(rounded[s] + hamming(x, s) for s in centres) for x in points}
    error = max(abs(f[x] - g[x]) for x in points)
    lipschitz = all(abs(g[x] - g[y]) <= hamming(x, y) for x in points for y in points)
    return {
        "cube_width": w,
        "cover_radius": radius,
        "cover_centres": len(centres),
        "advertised_error": delta,
        "observed_error": error,
        "extension_is_one_lipschitz": lipschitz,
        "extension_within_bound": error <= delta,
    }


def main() -> None:
    report = {
        "maxcut": maxcut_checks(),
        "pairwise_and_metric": pairwise_and_metric_checks(),
        "lipschitz_extension": lipschitz_check(),
    }
    report["all_checks_passed"] = all(
        (
            report["maxcut"]["global_flip_symmetry"],
            report["maxcut"]["pointwise_gluing"],
            report["maxcut"]["orbit_exposure"],
            report["pairwise_and_metric"]["selector_realization"],
            report["pairwise_and_metric"]["context_sup_metric"],
            report["pairwise_and_metric"]["shape_metric"],
            report["lipschitz_extension"]["extension_is_one_lipschitz"],
            report["lipschitz_extension"]["extension_within_bound"],
        )
    )
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
