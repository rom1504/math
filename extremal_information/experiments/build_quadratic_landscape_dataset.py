#!/usr/bin/env python3
"""Build an exact order-eight laboratory of quadratic sign landscapes.

Every signing can be switched so its first row is positive.  The remaining
negative edges form a graph on seven vertices, and NetworkX's graph atlas
contains one representative of every such graph-isomorphism class.  We group
the 1044 rooted gauges by their exact energy--energy--overlap histogram and
write one representative of each of the 243 resulting classes.

Integer energy, trace, near-cap, and one-vertex-response data are exact.  No
claim is made that the pair signature is a complete invariant at arbitrary
orders.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from itertools import product
from pathlib import Path

import networkx as nx
import numpy as np


def histogram(values: list[int]) -> list[list[int]]:
    return [[int(k), int(v)] for k, v in sorted(Counter(values).items())]


def matrix_from_graph(order: int, graph: nx.Graph) -> np.ndarray:
    matrix = np.ones((order, order), dtype=np.int64)
    np.fill_diagonal(matrix, 0)
    for i, j in graph.edges():
        matrix[i + 1, j + 1] = matrix[j + 1, i + 1] = -1
    return matrix


def pair_signature(
    energy: np.ndarray, overlap: np.ndarray, order: int
) -> bytes:
    edge_count = order * (order - 1) // 2
    energy_bins = 2 * edge_count + 1
    overlap_bins = 2 * order + 1
    index = (
        (energy[:, None] + edge_count) * energy_bins * overlap_bins
        + (energy[None, :] + edge_count) * overlap_bins
        + overlap
        + order
    ).ravel()
    return np.bincount(
        index, minlength=energy_bins * energy_bins * overlap_bins
    ).tobytes()


def record_for_graph(
    graph: nx.Graph,
    bucket_size: int,
    signature: bytes,
    spins: np.ndarray,
    overlap: np.ndarray,
) -> dict:
    order = spins.shape[1]
    matrix = matrix_from_graph(order, graph)
    energy = (
        np.einsum("bi,ij,bj->b", spins, matrix, spins, optimize=True) // 2
    ).astype(np.int64)
    absolute = np.abs(energy)
    cap = int(np.max(absolute))
    extension_caps = np.max(np.abs(energy[:, None] + overlap), axis=0)
    powers = {}
    matrix_power = np.eye(order, dtype=np.int64)
    for k in range(1, 9):
        matrix_power = matrix_power @ matrix
        powers[str(k)] = int(np.trace(matrix_power))
    return {
        "pair_signature_sha256": hashlib.sha256(signature).hexdigest(),
        "root_gauge_graph6": nx.to_graph6_bytes(
            graph, header=False
        ).decode().strip(),
        "root_gauge_graphs_in_signature_bucket": bucket_size,
        "matrix": matrix.astype(int).tolist(),
        "energy_histogram": histogram(energy.astype(int).tolist()),
        "cap": cap,
        "positive_maximum": int(np.max(energy)),
        "negative_minimum": int(np.min(energy)),
        "absolute_near_cap_counts": {
            str(gap): int(np.sum(absolute >= cap - gap)) for gap in (0, 2, 4)
        },
        "trace_powers_1_through_8": powers,
        "sorted_root_gauge_row_sums": sorted(
            matrix.sum(axis=1).astype(int).tolist()
        ),
        "one_vertex_cap_histogram": histogram(
            extension_caps.astype(int).tolist()
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--order", type=int, default=8, choices=[4, 5, 6, 7, 8])
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "extremal_information/experiments/quadratic_landscape_order8.json"
        ),
    )
    args = parser.parse_args()
    order = args.order
    spins = np.array(list(product((-1, 1), repeat=order)), dtype=np.int8)
    overlap = (spins @ spins.T).astype(np.int64)
    graphs = [g for g in nx.graph_atlas_g() if len(g) == order - 1]

    buckets: defaultdict[bytes, list[nx.Graph]] = defaultdict(list)
    for graph in graphs:
        matrix = matrix_from_graph(order, graph)
        energy = (
            np.einsum("bi,ij,bj->b", spins, matrix, spins, optimize=True) // 2
        ).astype(np.int64)
        buckets[pair_signature(energy, overlap, order)].append(graph)

    records = [
        record_for_graph(
            members[0], len(members), signature, spins, overlap
        )
        for signature, members in buckets.items()
    ]
    records.sort(key=lambda item: (item["cap"], item["pair_signature_sha256"]))
    minimum_cap = min(item["cap"] for item in records)
    for item in records:
        item["is_minimum_cap_class"] = item["cap"] == minimum_cap

    result = {
        "schema": "extremal-information-quadratic-landscape-dataset-v1",
        "order": order,
        "root_gauge_unlabeled_graphs": len(graphs),
        "pair_signature_classes": len(records),
        "minimum_cap": minimum_cap,
        "records": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {key: value for key, value in result.items() if key != "records"},
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
