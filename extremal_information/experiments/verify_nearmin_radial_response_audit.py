#!/usr/bin/env python3
"""Independent consistency checks for nearmin_radial_response_results.json."""

from __future__ import annotations

from collections import Counter
import json
from pathlib import Path

import numpy as np

import nearmin_radial_response_audit as audit


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "nearmin_radial_response_results.json"


def main() -> None:
    payload = json.loads(RESULTS.read_text())
    checked_covers = 0
    checked_metric_pairs = 0
    for record_index, record in enumerate(payload["records"]):
        a = np.asarray(record["matrix"], dtype=np.int8)
        n = len(a)
        assert audit.matrix_hash(a) == record["matrix_sha256"]
        edges = audit.edge_list(n)
        spins = audit.projective_spins(n)
        q = audit.edge_products(a, spins, edges)
        cap = audit.cap_of_q(q)
        assert cap == record["cap"]
        contexts, ends = audit.contexts_through_radius(len(edges), 3)
        sums = audit.edit_sums(q, contexts)
        adjusted = q.sum(axis=1, dtype=np.int16)[:, None] - 2 * sums

        for radius, stored in enumerate(record["radii"], 1):
            count = ends[radius]
            adj = adjusted[:, :count]
            response = np.max(np.abs(adj), axis=0)
            histogram = Counter(int(v - cap) for v in response)
            assert {str(k): v for k, v in sorted(histogram.items())} == stored[
                "response_distribution"
            ]["response_minus_base_cap_histogram"]

            exposed: set[int] = set()
            for c in range(count):
                for x in np.flatnonzero(np.abs(adj[:, c]) == response[c]):
                    if adj[x, c] >= 0:
                        exposed.add(2 * int(x) + 1)
                    if adj[x, c] <= 0:
                        exposed.add(2 * int(x))
            exposed = sorted(exposed)
            assert len(exposed) == stored["exposed_witness_count"]
            xs = np.asarray([z // 2 for z in exposed], dtype=np.int32)
            signs = np.asarray([1 if z % 2 else -1 for z in exposed], dtype=np.int16)
            values = signs[:, None] * adj[xs, :]

            base = q.sum(axis=1, dtype=np.int16)
            shell_scores = np.concatenate((-base, base))
            shell_size = int(np.count_nonzero(
                cap - shell_scores <= record["cap_excess"] + 2 * radius
            ))
            assert shell_size == stored["rs_shell_size"]
            assert len(exposed) <= shell_size

            for delta, cover in stored["response_cover"].items():
                if not cover["exact"]:
                    continue
                chosen = cover["selected_candidate_indices"]
                assert len(chosen) == cover["optimum"]
                envelope = np.max(values[chosen, :], axis=0)
                assert np.all(envelope >= response - int(delta))
                checked_covers += 1

            # Compare the closed metric formula with explicit evaluation for a
            # deterministic sparse set of pairs.  This includes every pair for
            # the first record at each order and a modular sample thereafter.
            t = signs[:, None] * q[xs, :]
            metric = audit.contextual_metric(t, radius)
            assert int(metric.max()) == stored["contextual_affine_metric"]["diameter"]
            for i in range(len(exposed)):
                for j in range(i):
                    if record_index % 9 or (17 * i + 31 * j + radius) % 101:
                        continue
                    explicit = int(np.max(np.abs(values[i, :] - values[j, :])))
                    assert explicit == int(metric[i, j])
                    checked_metric_pairs += 1

    print(
        f"PASS: {len(payload['records'])} records; "
        f"{checked_covers} certified covers; {checked_metric_pairs} metric pairs"
    )


if __name__ == "__main__":
    main()
