#!/usr/bin/env python3
"""Analyze all equal splits of a fixed parent into exact child minimizers."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np

from exact_mn_milp import projective_spins, stable_matrix_hash
from conference_double_construction import double_conference


def cap(matrix: np.ndarray, spins: np.ndarray) -> int:
    energies = np.einsum("bi,ij,bj->b", spins, matrix, spins) // 2
    return int(np.max(np.abs(energies)))


def rooted_gauge_key(matrix: np.ndarray) -> bytes:
    switches = np.ones(len(matrix), dtype=np.int8)
    switches[1:] = matrix[0, 1:]
    return (switches[:, None] * matrix * switches[None, :]).tobytes()


def signed_permutation_orbit(reference: np.ndarray, sign: int) -> set[bytes]:
    orbit: set[bytes] = set()
    for permutation in itertools.permutations(range(len(reference))):
        permuted = sign * reference[np.ix_(permutation, permutation)]
        orbit.add(rooted_gauge_key(permuted))
    return orbit


def add_orbit_class(
    classes: list[dict[str, set[bytes]]], reference: np.ndarray
) -> int:
    classes.append(
        {
            "plus": signed_permutation_orbit(reference, 1),
            "minus": signed_permutation_orbit(reference, -1),
        }
    )
    return len(classes) - 1


def classify(
    classes: list[dict[str, set[bytes]]], matrix: np.ndarray
) -> tuple[int, str]:
    key = rooted_gauge_key(matrix.astype(np.int8))
    for index, orbit in enumerate(classes):
        in_plus = key in orbit["plus"]
        in_minus = key in orbit["minus"]
        if in_plus or in_minus:
            sign = "+/-" if in_plus and in_minus else "+" if in_plus else "-"
            return index, sign
    index = add_orbit_class(classes, matrix.astype(np.int8))
    return index, "+"


def summary(values: list[float | int]) -> dict[str, float | int]:
    array = np.asarray(values, dtype=float)
    return {
        "minimum": float(array.min()),
        "maximum": float(array.max()),
        "mean": float(array.mean()),
        "distinct_count": int(len(np.unique(np.round(array, 12)))),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path)
    parser.add_argument("representative", type=Path)
    parser.add_argument("--representative-class", type=int)
    parser.add_argument("--matrix-key", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = json.loads(args.source.read_text())
    representative_payload = json.loads(args.representative.read_text())
    representative_data = (
        representative_payload["matrix"]
        if args.representative_class is None
        else representative_payload["classes"][args.representative_class][
            "representative_matrix"
        ]
    )
    representative = np.asarray(representative_data, dtype=np.int8)
    parent = (
        double_conference(representative)
        if args.matrix_key == "double_representative"
        else np.asarray(source[args.matrix_key], dtype=np.int8)
    )
    child_order = len(representative)
    if len(parent) != 2 * child_order:
        raise ValueError("the parent order is not twice the child order")
    spins = projective_spins(child_order).astype(np.int64)
    child_cap = cap(representative, spins)
    parent_spins = projective_spins(2 * child_order).astype(np.int64)
    parent_cap = cap(parent.astype(np.int64), parent_spins)
    classes: list[dict[str, set[bytes]]] = []
    add_orbit_class(classes, representative)

    correlations: list[float] = []
    active_counts: list[int] = []
    bridge_spectra: Counter[tuple[float, ...]] = Counter()
    equivalence_types: Counter[str] = Counter()
    class_spectrum_types: Counter[str] = Counter()
    exact_count = 0
    for subset in itertools.combinations(range(2 * child_order), child_order):
        if 0 not in subset:
            continue
        complement = tuple(
            vertex for vertex in range(2 * child_order) if vertex not in subset
        )
        left = parent[np.ix_(subset, subset)].astype(np.int64)
        right = parent[np.ix_(complement, complement)].astype(np.int64)
        if cap(left, spins) != child_cap or cap(right, spins) != child_cap:
            continue
        exact_count += 1
        bridge = parent[np.ix_(subset, complement)].astype(np.int64)
        left_energy = np.einsum("bi,ij,bj->b", spins, left, spins) // 2
        right_energy = np.einsum("bi,ij,bj->b", spins, right, spins) // 2
        internal = np.abs(left_energy[:, None] + right_energy[None, :])
        cross = np.abs(spins @ bridge @ spins.T)
        slack = parent_cap - internal - cross
        if int(slack.min()) < 0:
            raise AssertionError(int(slack.min()))
        correlations.append(float(np.corrcoef(internal.ravel(), cross.ravel())[0, 1]))
        active_counts.append(int(np.count_nonzero(slack == 0)))
        singular = tuple(np.round(np.linalg.svd(bridge, compute_uv=False), 10))
        bridge_spectra[singular] += 1

        left_class, left_sign = classify(classes, left)
        right_class, right_sign = classify(classes, right)
        equivalence_types[
            f"class{left_class}{left_sign},class{right_class}{right_sign}"
        ] += 1
        class_spectrum_types[
            f"class{left_class},class{right_class}|"
            + ",".join(f"{value:.10g}" for value in singular)
        ] += 1

    expected_from_source = source.get("exact_7_7_unordered_partition_count")
    if expected_from_source is not None and exact_count != expected_from_source:
        raise AssertionError((exact_count, expected_from_source))
    unordered_count = math.comb(2 * child_order, child_order) // 2
    output = {
        "schema": "quadratic-signing-equal-split-partition-analysis-v1",
        "classification": "exhaustive exact finite structural analysis",
        "source": str(args.source),
        "matrix_key": args.matrix_key,
        "representative": str(args.representative),
        "child_order": child_order,
        "child_cap": child_cap,
        "parent_order": 2 * child_order,
        "parent_cap": parent_cap,
        "parent_matrix_sha256": stable_matrix_hash(parent),
        "parent_matrix": [[int(value) for value in row] for row in parent],
        "unordered_partition_count": unordered_count,
        "exact_child_partition_count": exact_count,
        "signed_permutation_class_count": len(classes),
        "class_orbit_sizes": [
            {
                "plus": len(orbit["plus"]),
                "minus": len(orbit["minus"]),
                "intersection": len(orbit["plus"] & orbit["minus"]),
                "canonical_orbit_sha256": hashlib.sha256(
                    min(orbit["plus"] | orbit["minus"])
                ).hexdigest(),
            }
            for orbit in classes
        ],
        "block_equivalence_type_counts": dict(sorted(equivalence_types.items())),
        "class_and_bridge_spectrum_counts": dict(
            sorted(class_spectrum_types.items())
        ),
        "margin_correlation": summary(correlations),
        "active_constraint_count": summary(active_counts),
        "bridge_singular_spectrum_counts": {
            ",".join(f"{value:.10g}" for value in spectrum): count
            for spectrum, count in sorted(bridge_spectra.items())
        },
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(
        f"exact={exact_count}/{unordered_count} equivalence={dict(equivalence_types)} "
        f"spectral_types={len(bridge_spectra)}"
    )
    print(
        f"correlation={output['margin_correlation']} "
        f"active={output['active_constraint_count']}"
    )
    print(f"wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
