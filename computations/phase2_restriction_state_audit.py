#!/usr/bin/env python3
"""Test bounded principal-restriction states on Paley subset orbits.

For every saved child-orbit representative, compute the complete histogram of
switching/permutation/global-negation classes of all t-vertex principal
restrictions.  Combine these profiles cumulatively with the child's complete
spectrum and ask whether the state determines its Boolean cap.
"""

from __future__ import annotations

import argparse
import itertools
import json
from collections import Counter
from pathlib import Path

import numpy as np

from conference_prime_square import PrimeSquare, stable_hash


class DisjointSet:
    def __init__(self, n: int):
        self.parent = list(range(n))

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: int, b: int) -> None:
        a, b = self.find(a), self.find(b)
        if a != b:
            self.parent[max(a, b)] = min(a, b)


def signing_from_code(code: int, n: int) -> np.ndarray:
    matrix = np.ones((n, n), dtype=np.int8)
    np.fill_diagonal(matrix, 0)
    bit = 0
    for i in range(1, n):
        for j in range(i + 1, n):
            if code & (1 << bit):
                matrix[i, j] = matrix[j, i] = -1
            bit += 1
    return matrix


def root_gauge_code(matrix: np.ndarray) -> int:
    switches = np.ones(len(matrix), dtype=np.int8)
    switches[1:] = matrix[0, 1:]
    gauged = switches[:, None] * matrix * switches[None, :]
    code = 0
    bit = 0
    for i in range(1, len(matrix)):
        for j in range(i + 1, len(matrix)):
            if gauged[i, j] == -1:
                code |= 1 << bit
            bit += 1
    return code


def class_map(n: int) -> tuple[list[int], int]:
    bits = (n - 1) * (n - 2) // 2
    count = 1 << bits
    dsu = DisjointSet(count)
    adjacent = []
    for i in range(n - 1):
        permutation = list(range(n))
        permutation[i], permutation[i + 1] = permutation[i + 1], permutation[i]
        adjacent.append(permutation)
    full_mask = count - 1
    for code in range(count):
        matrix = signing_from_code(code, n)
        dsu.union(code, code ^ full_mask)
        for permutation in adjacent:
            moved = matrix[np.ix_(permutation, permutation)]
            dsu.union(code, root_gauge_code(moved))
    roots = [dsu.find(code) for code in range(count)]
    labels = {root: label for label, root in enumerate(sorted(set(roots)))}
    return [labels[root] for root in roots], len(labels)


def profile(matrix: np.ndarray, size: int, labels: list[int]) -> tuple[int, ...]:
    counts: Counter[int] = Counter()
    for vertices in itertools.combinations(range(len(matrix)), size):
        child = matrix[np.ix_(vertices, vertices)]
        counts[labels[root_gauge_code(child)]] += 1
    return tuple(counts[index] for index in range(max(labels) + 1))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("orbit_result", type=Path)
    parser.add_argument("--p", type=int, required=True)
    parser.add_argument("--max-restriction-size", type=int, default=7)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.orbit_result.read_text())
    field = PrimeSquare(args.p)
    matrix = field.conference()
    if stable_hash(matrix) != source["conference_matrix_sha256"]:
        raise AssertionError("conference matrix hash mismatch")
    records = source["orbit_spectral_records"]

    maps = {}
    class_counts = {}
    for size in range(4, args.max_restriction_size + 1):
        maps[size], class_counts[size] = class_map(size)
        print(f"size={size} classes={class_counts[size]}", flush=True)

    enriched = []
    for index, record in enumerate(records, start=1):
        vertices = tuple(record["vertices"])
        child = matrix[np.ix_(vertices, vertices)]
        profiles = {
            str(size): profile(child, size, maps[size]) for size in maps
        }
        enriched.append({**record, "restriction_class_histograms": profiles})
        if index % 50 == 0:
            print(f"children={index}/{len(records)}", flush=True)

    audits = {}
    restriction_only_audits = {}
    for through in range(4, args.max_restriction_size + 1):
        groups: dict[tuple[object, ...], list[dict[str, object]]] = {}
        for record in enriched:
            signature = (
                tuple(record["characteristic_coefficients"]),
                *(tuple(record["restriction_class_histograms"][str(size)])
                  for size in range(4, through + 1)),
            )
            groups.setdefault(signature, []).append(record)
        ambiguous = [items for items in groups.values()
                     if len({int(item["cap"]) for item in items}) > 1]
        examples = []
        for items in ambiguous[:5]:
            by_cap = {}
            for item in items:
                by_cap.setdefault(int(item["cap"]), item)
            examples.append([
                {"cap": cap, "vertices": item["vertices"],
                 "orbit_size": item["orbit_size"]}
                for cap, item in sorted(by_cap.items())
            ])
        audits[str(through)] = {
            "state": f"full spectrum plus complete restriction-class histograms for sizes 4..{through}",
            "distinct_states": len(groups),
            "ambiguous_states": len(ambiguous),
            "ambiguous_cap_sets": dict(Counter(
                ",".join(map(str, sorted({int(item["cap"]) for item in items})))
                for items in ambiguous)),
            "cap_determined": not ambiguous,
            "ambiguous_examples": examples,
        }

        restriction_groups: dict[tuple[object, ...], list[dict[str, object]]] = {}
        for record in enriched:
            signature = tuple(
                tuple(record["restriction_class_histograms"][str(size)])
                for size in range(4, through + 1)
            )
            restriction_groups.setdefault(signature, []).append(record)
        restriction_ambiguous = [
            items for items in restriction_groups.values()
            if len({int(item["cap"]) for item in items}) > 1]
        restriction_only_audits[str(through)] = {
            "state": f"complete restriction-class histograms for sizes 4..{through}, without spectrum",
            "distinct_states": len(restriction_groups),
            "ambiguous_states": len(restriction_ambiguous),
            "ambiguous_cap_sets": dict(Counter(
                ",".join(map(str, sorted({int(item["cap"]) for item in items})))
                for items in restriction_ambiguous)),
            "cap_determined": not restriction_ambiguous,
        }

    output = {
        "schema": "quadratic-signing-phase2-restriction-state-audit-v1",
        "classification": (
            "exhaustive invariant audit over every certified Paley subset "
            "orbit; switching/permutation/global-negation classes exact"
        ),
        "source": str(args.orbit_result),
        "p": args.p,
        "conference_order": len(matrix),
        "child_order": len(records[0]["vertices"]),
        "child_orbits": len(records),
        "restriction_class_counts": {
            str(size): count for size, count in class_counts.items()},
        "cumulative_state_audits": audits,
        "restriction_only_state_audits": restriction_only_audits,
        "enriched_orbit_records": enriched,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
