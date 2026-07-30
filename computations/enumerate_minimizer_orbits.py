#!/usr/bin/env python3
"""Exhaustively classify small-order minimizers up to basic equivalence."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from collections import Counter
from pathlib import Path

import numpy as np

from analyze_equal_split_partitions import add_orbit_class, classify
from exact_mn_milp import exact_profile, projective_spins, stable_matrix_hash


def orbit_id(orbit: dict[str, set[bytes]]) -> str:
    return hashlib.sha256(min(orbit["plus"] | orbit["minus"])).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("representative", type=Path)
    parser.add_argument("--partition-analysis", type=Path)
    parser.add_argument("--require-cover", action="store_true")
    parser.add_argument("--count-only", action="store_true")
    parser.add_argument("--batch-size", type=int, default=8192)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    representative_payload = json.loads(args.representative.read_text())
    representative = np.asarray(representative_payload["matrix"], dtype=np.int8)
    n = len(representative)
    target_cap = int(representative_payload["profile"]["M"])
    spins = projective_spins(n).astype(np.int16)
    internal_edges = tuple(itertools.combinations(range(1, n), 2))
    classes: list[dict[str, set[bytes]]] = []
    add_orbit_class(classes, representative)
    class_representatives = [representative.copy()]
    class_counts: Counter[int] = Counter()
    minimizer_count = 0
    minimizing_masks: list[int] = []
    edge_products = np.asarray(
        [spins[:, i] * spins[:, j] for i, j in internal_edges], dtype=np.int16
    ).T
    root_energy = np.sum(spins[:, 1:], axis=1, dtype=np.int16)
    total_signings = 1 << len(internal_edges)
    bit_positions = np.arange(len(internal_edges), dtype=np.uint64)
    for start in range(0, total_signings, args.batch_size):
        masks = np.arange(
            start, min(start + args.batch_size, total_signings), dtype=np.uint64
        )
        signs = 1 - 2 * ((masks[:, None] >> bit_positions) & 1).astype(np.int16)
        energies = signs @ edge_products.T + root_energy[None, :]
        feasible = np.max(np.abs(energies), axis=1) == target_cap
        minimizing_masks.extend(int(mask) for mask in masks[feasible])
    minimizer_count = len(minimizing_masks)
    for mask in minimizing_masks:
        matrix = np.zeros((n, n), dtype=np.int8)
        matrix[0, 1:] = matrix[1:, 0] = 1
        for bit, (i, j) in enumerate(internal_edges):
            matrix[i, j] = matrix[j, i] = -1 if (mask >> bit) & 1 else 1
        if args.count_only:
            continue
        previous_class_count = len(classes)
        class_index, _ = classify(classes, matrix)
        if len(classes) > previous_class_count:
            class_representatives.append(matrix.copy())
        class_counts[class_index] += 1

    class_rows = []
    for index, orbit in enumerate(classes if not args.count_only else []):
        union = orbit["plus"] | orbit["minus"]
        class_rows.append(
            {
                "class": index,
                "canonical_orbit_sha256": orbit_id(orbit),
                "root_gauged_labeled_count": class_counts[index],
                "enumerated_orbit_union_size": len(union),
                "self_complementary": bool(orbit["plus"] & orbit["minus"]),
                "representative_matrix_sha256": stable_matrix_hash(
                    class_representatives[index]
                ),
                "representative_matrix": [
                    [int(value) for value in row]
                    for row in class_representatives[index]
                ],
                "representative_profile": exact_profile(
                    class_representatives[index]
                ),
            }
        )
        if class_counts[index] != len(union):
            raise AssertionError((index, class_counts[index], len(union)))
    conference_ids: list[str] | None = None
    covers_conference_classes: bool | None = None
    if args.partition_analysis and not args.count_only:
        partition_payload = json.loads(args.partition_analysis.read_text())
        conference_ids = sorted(
            row["canonical_orbit_sha256"]
            for row in partition_payload["class_orbit_sizes"]
        )
        covers_conference_classes = conference_ids == sorted(
            row["canonical_orbit_sha256"] for row in class_rows
        )
        if args.require_cover and not covers_conference_classes:
            raise AssertionError("conference and exhaustive class sets differ")
    output = {
        "schema": "quadratic-signing-minimizer-orbits-v1",
        "classification": "exhaustive exact finite enumeration",
        "normalization": f"root switching gauge fixes all {n-1} root edges to +1",
        "order": n,
        "target_cap": target_cap,
        "root_gauged_signing_count": 1 << len(internal_edges),
        "minimizing_signing_count": minimizer_count,
        "signed_permutation_and_global_sign_class_count": (
            None if args.count_only else len(classes)
        ),
        "classes": class_rows,
        "conference_class_ids": conference_ids,
        "conference_splits_cover_all_minimizer_classes": covers_conference_classes,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"n={n} enumerated={output['root_gauged_signing_count']} "
        f"cap{target_cap}={minimizer_count} "
        f"classes={'not-enumerated' if args.count_only else len(classes)} "
        f"counts={dict(class_counts)}"
    )
    print(f"conference covers classes={covers_conference_classes}")
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
