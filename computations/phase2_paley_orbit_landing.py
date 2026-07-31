#!/usr/bin/env python3
"""Orbit-complete balanced landing audit for ``PC(p^2+1)``.

The projective semilinear group acts on the Paley two-graph.  Its action sends
every principal signing to a switching-equivalent signing, up to global
negation, so one exact cap evaluation per subset orbit suffices.  This script
constructs that action directly over GF(p^2), enumerates every balanced-subset
orbit, and checks one representative with the batch Gray evaluator.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import subprocess
import time
from collections import Counter, deque
from pathlib import Path

import numpy as np

from conference_prime_square import PrimeSquare, stable_hash


def compose(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    """Return left after right."""
    return tuple(left[right[i]] for i in range(len(left)))


def group_closure(generators: list[tuple[int, ...]]) -> list[tuple[int, ...]]:
    identity = tuple(range(len(generators[0])))
    seen = {identity}
    queue = deque([identity])
    ordered = []
    while queue:
        item = queue.popleft()
        ordered.append(item)
        for generator in generators:
            nxt = compose(generator, item)
            if nxt not in seen:
                seen.add(nxt)
                queue.append(nxt)
    return ordered


def evaluate_batch(executable: Path, matrix: np.ndarray,
                   subsets: list[tuple[int, ...]]) -> list[tuple[int, int, int]]:
    n, k = len(matrix), len(subsets[0])
    chunks = [f"{n} {k} {len(subsets)}\n"]
    chunks.append(" ".join(map(str, matrix.reshape(-1))) + "\n")
    chunks.extend(" ".join(map(str, subset)) + "\n" for subset in subsets)
    result = subprocess.run([str(executable)], input="".join(chunks), text=True,
                            capture_output=True, check=True)
    rows = [tuple(map(int, line.split())) for line in result.stdout.splitlines()]
    if len(rows) != len(subsets):
        raise AssertionError((len(rows), len(subsets), result.stderr))
    return rows


def exact_spectral_data(matrix: np.ndarray) -> tuple[list[int], list[int]]:
    """Return exact power traces and monic characteristic coefficients."""
    n = len(matrix)
    integer = matrix.astype(np.int64)
    power = np.eye(n, dtype=np.int64)
    traces = []
    for _ in range(1, n + 1):
        power = power @ integer
        traces.append(int(np.trace(power)))
    coefficients = [1]
    for k in range(1, n + 1):
        numerator = sum(coefficients[k - i] * traces[i - 1]
                        for i in range(1, k + 1))
        if numerator % k:
            raise AssertionError((k, numerator))
        coefficients.append(-numerator // k)
    # Exact Cayley--Hamilton check protects against integer overflow and a
    # faulty Newton-identity implementation.
    accumulator = np.zeros((n, n), dtype=object)
    integer_object = matrix.astype(object)
    identity_object = np.eye(n, dtype=object)
    for coefficient in coefficients:
        accumulator = accumulator @ integer_object + coefficient * identity_object
    if np.any(accumulator):
        raise AssertionError("Cayley--Hamilton check failed")
    return traces, coefficients


def canonical_signed_prefix(traces: list[int], degree: int) -> tuple[int, ...]:
    direct = tuple(traces[i - 1] for i in range(3, degree + 1))
    negated = tuple(((-1) ** i) * traces[i - 1] for i in range(3, degree + 1))
    return min(direct, negated)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--evaluator", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    field = PrimeSquare(args.p)
    q = args.p ** 2
    n = q + 1
    if n % 2:
        raise ValueError("conference order must be even")
    k = n // 2
    elements = list(field.elements)
    index = {value: i + 1 for i, value in enumerate(elements)}
    zero, one, t = (0, 0), (1, 0), (0, 1)

    def add(x: tuple[int, int], y: tuple[int, int]) -> tuple[int, int]:
        return ((x[0] + y[0]) % args.p, (x[1] + y[1]) % args.p)

    def inverse(x: tuple[int, int]) -> tuple[int, int]:
        if x == zero:
            raise ZeroDivisionError
        return field.power(x, q - 2)

    # Index 0 denotes infinity.  Translations by a GF(p)-basis, one arbitrary
    # nonzero scaling, inversion, and Frobenius generate PGammaL(2,q).
    generators: list[tuple[int, ...]] = []
    for shift in (one, t):
        generators.append((0,) + tuple(index[add(x, shift)] for x in elements))
    nonsquare = next(x for x in elements if field.character(x) == -1)
    generators.append((0,) + tuple(index[field.multiply(nonsquare, x)] for x in elements))
    inversion = [index[zero], 0]
    inversion.extend(index[inverse(x)] for x in elements[1:])
    generators.append(tuple(inversion))
    generators.append((0,) + tuple(index[field.power(x, args.p)] for x in elements))

    matrix = field.conference()

    # Independently verify on the actual integer signing that every generator
    # acts by switching and, when needed, global negation.  Generator checks
    # imply cap invariance for the generated group.
    generator_actions = []
    for permutation in generators:
        moved = matrix[np.ix_(permutation, permutation)].astype(np.int64)
        certificate = None
        for global_sign in (1, -1):
            switches = np.ones(n, dtype=np.int64)
            switches[1:] = moved[0, 1:] * matrix[0, 1:] * global_sign
            reconstructed = global_sign * switches[:, None] * matrix * switches[None, :]
            if np.array_equal(moved, reconstructed):
                certificate = {
                    "global_sign": global_sign,
                    "negative_switch_vertices": [
                        i for i, value in enumerate(switches) if value == -1],
                }
                break
        if certificate is None:
            raise AssertionError("generator is not switching/global-negation invariant")
        generator_actions.append(certificate)

    started = time.monotonic()
    group = group_closure(generators)
    expected_group_order = 2 * q * (q * q - 1)
    if len(group) != expected_group_order:
        raise AssertionError((len(group), expected_group_order))
    print(f"group order {len(group)}", flush=True)

    # A byte per subset mask keeps the exhaustive orbit marking simple and
    # still uses only 2^(q+1) bytes at the intended order 26 experiment.
    if n > 27:
        raise ValueError("explicit orbit bitmap is restricted to order <=27")
    visited = bytearray(1 << n)
    representatives: list[tuple[int, ...]] = []
    orbit_sizes: list[int] = []
    total_subsets = math.comb(n, k)
    marked = 0
    for subset in itertools.combinations(range(n), k):
        mask = sum(1 << vertex for vertex in subset)
        if visited[mask]:
            continue
        orbit_masks = set()
        for permutation in group:
            orbit_masks.add(sum(1 << permutation[vertex] for vertex in subset))
        for image in orbit_masks:
            if not visited[image]:
                visited[image] = 1
                marked += 1
        representatives.append(subset)
        orbit_sizes.append(len(orbit_masks))
        if len(representatives) % 50 == 0:
            print(f"orbits={len(representatives)} marked={marked}/{total_subsets}",
                  flush=True)
    if marked != total_subsets:
        raise AssertionError((marked, total_subsets))

    profiles = evaluate_batch(args.evaluator.resolve(), matrix, representatives)
    caps = [row[0] for row in profiles]
    spectral_records = []
    for subset, cap, orbit_size in zip(representatives, caps, orbit_sizes):
        child = matrix[np.ix_(subset, subset)]
        traces, coefficients = exact_spectral_data(child)
        negative_coefficients = [
            value * ((-1) ** degree)
            for degree, value in enumerate(coefficients)
        ]
        spectral_records.append({
            "vertices": list(subset),
            "orbit_size": orbit_size,
            "cap": cap,
            "power_traces_3_through_n": traces[2:],
            "characteristic_coefficients": list(
                min(tuple(coefficients), tuple(negative_coefficients))),
        })

    spectral_state_audit = {}
    for degree in range(3, k + 1):
        states: dict[tuple[int, ...], set[int]] = {}
        for record in spectral_records:
            child_traces = [0, k * (k - 1),
                            *record["power_traces_3_through_n"]]
            signature = canonical_signed_prefix(child_traces, degree)
            states.setdefault(signature, set()).add(int(record["cap"]))
        ambiguous = [values for values in states.values() if len(values) > 1]
        spectral_state_audit[str(degree)] = {
            "state": f"power traces through degree {degree}, modulo global negation",
            "distinct_states": len(states),
            "ambiguous_states": len(ambiguous),
            "ambiguous_cap_sets": dict(Counter(
                ",".join(map(str, sorted(values))) for values in ambiguous)),
            "cap_determined": not ambiguous,
        }
    characteristic_states: dict[tuple[int, ...], list[dict[str, object]]] = {}
    for record in spectral_records:
        signature = tuple(record["characteristic_coefficients"])
        characteristic_states.setdefault(signature, []).append(record)
    ambiguous_characteristic = [
        records for records in characteristic_states.values()
        if len({int(record["cap"]) for record in records}) > 1]
    ambiguous_characteristic_examples = []
    for records in ambiguous_characteristic[:5]:
        by_cap: dict[int, dict[str, object]] = {}
        for record in records:
            by_cap.setdefault(int(record["cap"]), record)
        ambiguous_characteristic_examples.append({
            "characteristic_coefficients_modulo_global_negation":
                records[0]["characteristic_coefficients"],
            "representatives": [
                {"vertices": item["vertices"], "orbit_size": item["orbit_size"],
                 "cap": item["cap"]}
                for _, item in sorted(by_cap.items())
            ],
        })
    spectral_state_audit["full_characteristic_polynomial"] = {
        "state": "complete spectrum modulo global negation",
        "distinct_states": len(characteristic_states),
        "ambiguous_states": len(ambiguous_characteristic),
        "ambiguous_cap_sets": dict(Counter(
            ",".join(map(str, sorted({int(record["cap"]) for record in records})))
            for records in ambiguous_characteristic)),
        "ambiguous_examples": ambiguous_characteristic_examples,
        "cap_determined": not ambiguous_characteristic,
    }
    weighted_histogram: Counter[int] = Counter()
    orbit_histogram: Counter[int] = Counter()
    for cap, orbit_size in zip(caps, orbit_sizes):
        weighted_histogram[cap] += orbit_size
        orbit_histogram[cap] += 1
    best = min(caps)
    best_records = [
        {"vertices": list(representatives[i]), "orbit_size": orbit_sizes[i],
         "cap": profiles[i][0], "min_energy": profiles[i][1],
         "max_energy": profiles[i][2]}
        for i in range(len(caps)) if caps[i] == best
    ]
    serialized = json.dumps(
        [[list(s), orbit_sizes[i], list(profiles[i])]
         for i, s in enumerate(representatives)], separators=(",", ":"))
    output = {
        "schema": "quadratic-signing-phase2-paley-orbit-landing-v1",
        "classification": (
            "exhaustive balanced-subset orbit classification under the exact "
            "PGammaL action; exact Gray cap for every orbit representative"
        ),
        "field": f"GF({q})=GF({args.p})[t]/(t^2-{field.d})",
        "conference_order": n,
        "child_order": k,
        "conference_matrix_sha256": stable_hash(matrix),
        "group": "PGammaL(2,q)",
        "group_order": len(group),
        "group_order_formula": "2*q*(q^2-1)",
        "generator_action_certificates": generator_actions,
        "cap_invariance_reason": (
            "a fractional-linear map multiplies chi(x-y) by one global "
            "character and two vertex characters; Frobenius preserves chi"
        ),
        "total_balanced_subsets": total_subsets,
        "subset_orbits": len(representatives),
        "orbit_size_histogram": {
            str(size): count for size, count in sorted(Counter(orbit_sizes).items())},
        "minimum_cap": best,
        "minimum_cap_orbits": best_records,
        "cap_orbit_histogram": {
            str(cap): count for cap, count in sorted(orbit_histogram.items())},
        "cap_subset_histogram": {
            str(cap): count for cap, count in sorted(weighted_histogram.items())},
        "spectral_state_audit": spectral_state_audit,
        "orbit_spectral_records": spectral_records,
        "all_orbit_records_sha256": hashlib.sha256(serialized.encode()).hexdigest(),
        "evaluator_source": "computations/phase2_subset_caps_gray.cpp",
        "elapsed_seconds": time.monotonic() - started,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(f"orbits={len(representatives)} minimum_cap={best} wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
