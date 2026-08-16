#!/usr/bin/env python3
"""Exact finite tests for extremal summaries of Boolean landscapes.

The script verifies three claims used in the initial theory audit:

1. two four-bit codes have the same complete pair-distance data but different
   covering radii, and this separation tensorizes;
2. two order-eight quadratic signings have identical energy histograms but a
   different response to the same one-vertex coupling; and
3. through order eight, the exact joint histogram of
   (energy(x), energy(y), overlap(x,y)) has no collision with different
   one-vertex response multisets.

All enumeration is exact.  The graph-atlas census is complete because after
switching the first row of a signing to +1, its remaining negative edges form
an arbitrary graph on n-1 vertices.
"""

from __future__ import annotations

import argparse
import json
from collections import Counter, defaultdict
from itertools import product
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import networkx as nx
import numpy as np


CODE_C = (0b0000, 0b0011, 0b0101, 0b0110)
CODE_D = (0b0000, 0b0011, 0b0101, 0b1001)
QUADRATIC_MASK_A = 1466915
QUADRATIC_MASK_B = 1068688


def hamming(a: int, b: int) -> int:
    return bin(a ^ b).count("1")


def distance_enumerator(code: Sequence[int]) -> Dict[int, int]:
    return dict(
        sorted(
            Counter(hamming(a, b) for a in code for b in code).items()
        )
    )


def covering_radius(code: Sequence[int], dimension: int) -> int:
    return max(
        min(hamming(x, c) for c in code) for x in range(1 << dimension)
    )


def tensor_distance_enumerator(
    base: Dict[int, int], tensor_power: int
) -> Dict[int, int]:
    answer = {0: 1}
    for _ in range(tensor_power):
        updated: Counter[int] = Counter()
        for d1, c1 in answer.items():
            for d2, c2 in base.items():
                updated[d1 + d2] += c1 * c2
        answer = dict(updated)
    return dict(sorted(answer.items()))


def spin_cube(order: int) -> np.ndarray:
    return np.array(list(product((-1, 1), repeat=order)), dtype=np.int8)


def gauge_matrix_from_mask(order: int, mask: int) -> np.ndarray:
    matrix = np.ones((order, order), dtype=np.int8)
    np.fill_diagonal(matrix, 0)
    free_edges = [
        (i, j) for i in range(1, order) for j in range(i + 1, order)
    ]
    for bit, (i, j) in enumerate(free_edges):
        matrix[i, j] = matrix[j, i] = 1 if (mask >> bit) & 1 else -1
    return matrix


def gauge_matrix_from_graph(order: int, graph: nx.Graph) -> np.ndarray:
    matrix = np.ones((order, order), dtype=np.int8)
    np.fill_diagonal(matrix, 0)
    for i, j in graph.edges():
        matrix[i + 1, j + 1] = matrix[j + 1, i + 1] = -1
    return matrix


def energies(matrix: np.ndarray, spins: np.ndarray) -> np.ndarray:
    return (
        np.einsum("bi,ij,bj->b", spins, matrix, spins, optimize=True) // 2
    ).astype(np.int64)


def histogram(values: Iterable[int]) -> List[Tuple[int, int]]:
    return sorted((int(k), int(v)) for k, v in Counter(values).items())


def one_vertex_responses(
    energy: np.ndarray, spins: np.ndarray, overlaps: np.ndarray
) -> List[Tuple[int, int]]:
    # Columns range over every incident sign vector b.  Since the new spin can
    # change sign, max_x |H(x)+<b,x>| is the exact extended cap.
    caps = np.max(np.abs(energy[:, None] + overlaps), axis=0)
    return histogram(caps.tolist())


def pair_signature_bytes(
    energy: np.ndarray, overlaps: np.ndarray, order: int
) -> bytes:
    edge_count = order * (order - 1) // 2
    energy_bins = 2 * edge_count + 1
    overlap_bins = 2 * order + 1
    indices = (
        (energy[:, None] + edge_count) * energy_bins * overlap_bins
        + (energy[None, :] + edge_count) * overlap_bins
        + overlaps
        + order
    ).ravel()
    return np.bincount(
        indices, minlength=energy_bins * energy_bins * overlap_bins
    ).tobytes()


def code_collision(max_tensor_power: int) -> dict:
    enum_c = distance_enumerator(CODE_C)
    enum_d = distance_enumerator(CODE_D)
    assert enum_c == enum_d
    radius_c = covering_radius(CODE_C, 4)
    radius_d = covering_radius(CODE_D, 4)
    assert (radius_c, radius_d) == (2, 3)
    powers = []
    for r in range(1, max_tensor_power + 1):
        tensor_enum_c = tensor_distance_enumerator(enum_c, r)
        tensor_enum_d = tensor_distance_enumerator(enum_d, r)
        assert tensor_enum_c == tensor_enum_d
        powers.append(
            {
                "tensor_power": r,
                "dimension": 4 * r,
                "covering_radius_C": radius_c * r,
                "covering_radius_D": radius_d * r,
                "normalized_gap": (radius_d - radius_c) / 4,
                "distance_enumerator": tensor_enum_c,
            }
        )
    return {
        "C": list(CODE_C),
        "D": list(CODE_D),
        "ordered_distance_enumerator": enum_c,
        "covering_radius_C": radius_c,
        "covering_radius_D": radius_d,
        "tensor_powers": powers,
    }


def scalar_entropy_collision() -> dict:
    order = 8
    spins = spin_cube(order)
    overlaps = (spins @ spins.T).astype(np.int64)
    records = []
    for mask in (QUADRATIC_MASK_A, QUADRATIC_MASK_B):
        matrix = gauge_matrix_from_mask(order, mask)
        energy = energies(matrix, spins)
        b = -np.ones(order, dtype=np.int8)
        fixed_cap = int(np.max(np.abs(energy + spins @ b)))
        records.append(
            {
                "mask": mask,
                "matrix": matrix.astype(int).tolist(),
                "energy_histogram": histogram(energy.tolist()),
                "base_cap": int(np.max(np.abs(energy))),
                "all_minus_one_vertex_cap": fixed_cap,
                "one_vertex_response_histogram": one_vertex_responses(
                    energy, spins, overlaps
                ),
            }
        )
    assert records[0]["energy_histogram"] == records[1]["energy_histogram"]
    assert records[0]["base_cap"] == records[1]["base_cap"] == 14
    assert (
        records[0]["all_minus_one_vertex_cap"],
        records[1]["all_minus_one_vertex_cap"],
    ) == (16, 20)
    return {"order": order, "records": records}


def pair_overlap_census(min_order: int, max_order: int) -> list:
    census = []
    atlas = nx.graph_atlas_g()
    for order in range(min_order, max_order + 1):
        spins = spin_cube(order)
        overlaps = (spins @ spins.T).astype(np.int64)
        graphs = [g for g in atlas if len(g) == order - 1]
        buckets: defaultdict[bytes, list] = defaultdict(list)
        for atlas_index, graph in enumerate(graphs):
            energy = energies(gauge_matrix_from_graph(order, graph), spins)
            signature = pair_signature_bytes(energy, overlaps, order)
            response = one_vertex_responses(energy, spins, overlaps)
            buckets[signature].append((atlas_index, response))
        separating = []
        for representatives in buckets.values():
            response_types = {
                tuple(tuple(v) for v in response): atlas_index
                for atlas_index, response in representatives
            }
            if len(response_types) > 1:
                separating.append(sorted(response_types.values()))
        assert not separating
        census.append(
            {
                "order": order,
                "root_gauge_unlabeled_graphs": len(graphs),
                "exact_pair_energy_overlap_signatures": len(buckets),
                "signature_buckets_with_multiple_root_gauges": sum(
                    len(v) > 1 for v in buckets.values()
                ),
                "different_response_collisions": separating,
            }
        )
    return census


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "extremal_information/experiments/entropy_overlap_results.json"
        ),
    )
    parser.add_argument("--max-tensor-power", type=int, default=5)
    parser.add_argument("--min-census-order", type=int, default=4)
    parser.add_argument("--max-census-order", type=int, default=8)
    args = parser.parse_args()

    result = {
        "schema": "extremal-information-entropy-overlap-v1",
        "code_pair_overlap_collision": code_collision(args.max_tensor_power),
        "quadratic_scalar_entropy_collision": scalar_entropy_collision(),
        "quadratic_pair_overlap_census": pair_overlap_census(
            args.min_census_order, args.max_census_order
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
