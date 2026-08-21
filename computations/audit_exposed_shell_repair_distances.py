#!/usr/bin/env python3
"""Exact sparse-repair distances on exposed cap shells.

For a root-gauged signing, flipping a physical edge induces a generator of
the switching quotient.  Multi-source breadth-first search in this Cayley
graph therefore gives the minimum number of physical edge flips, modulo
switching, needed to enter a lower cap layer.  Permutation invariance lets us
fix the restricted vertex set to ``range(m)`` when conditioning on a bad
restriction incidence.

The default audit is exhaustive at order eight.  It extends the one-edge
audit in ``audit_exposed_shell_edge_repairs.py`` to the complete distance
distribution, both for arbitrary parent edges and for edges internal to the
restricted set.
"""

from __future__ import annotations

import argparse
import itertools
import json
import math
from pathlib import Path

import numpy as np

from audit_canonical_disorder_root_gauge import projected_masks, root_caps
from audit_exposed_shell_edge_repairs import edge_toggles


UNREACHABLE = np.iinfo(np.uint8).max


def validate_physical_generators(n: int) -> None:
    """Check the quotient generator formula by explicit flip and regauging."""
    expected = edge_toggles(n, range(n))
    observed: list[int] = []
    for flipped_i, flipped_j in itertools.combinations(range(n), 2):
        # Begin at the all-positive root gauge, flip one physical edge, and
        # switch every non-root vertex whose new root edge is negative.
        matrix = np.ones((n, n), dtype=np.int8)
        np.fill_diagonal(matrix, 0)
        matrix[flipped_i, flipped_j] *= -1
        matrix[flipped_j, flipped_i] *= -1
        switch = np.ones(n, dtype=np.int8)
        switch[1:] = matrix[0, 1:]
        matrix *= switch[:, None] * switch[None, :]
        assert np.all(matrix[0, 1:] == 1)
        mask = 0
        for bit, (i, j) in enumerate(itertools.combinations(range(1, n), 2)):
            if matrix[i, j] == -1:
                mask |= 1 << bit
        observed.append(mask)
    assert observed == expected, (observed, expected)


def quotient_distances(
    parent_caps: np.ndarray,
    target_cap: int,
    generators: list[int],
) -> np.ndarray:
    """Return exact Cayley-graph distance to ``cap <= target_cap``.

    An entry is ``UNREACHABLE`` precisely when its connected component under
    the allowed physical-edge generators contains no target signing.
    """
    state_count = len(parent_caps)
    if state_count > 1 << 32:
        raise ValueError("audit currently requires at most 2^32 quotient states")
    indices = np.arange(state_count, dtype=np.uint32)
    visited = parent_caps <= target_cap
    distance = np.full(state_count, UNREACHABLE, dtype=np.uint8)
    distance[visited] = 0
    depth = 0
    while True:
        reachable = visited.copy()
        for generator in generators:
            reachable |= visited[indices ^ np.uint32(generator)]
        new = reachable & ~visited
        if not np.any(new):
            return distance
        depth += 1
        if depth >= UNREACHABLE:
            raise RuntimeError("distance overflow")
        distance[new] = depth
        visited = reachable


def distance_histogram(distance: np.ndarray) -> dict[str, int]:
    values, counts = np.unique(distance, return_counts=True)
    return {
        ("unreachable" if int(value) == UNREACHABLE else str(int(value))): int(count)
        for value, count in zip(values, counts)
    }


def finite_radius(distance: np.ndarray) -> int | None:
    finite = distance[distance != UNREACHABLE]
    return int(np.max(finite)) if len(finite) else None


def audit(n: int, child_orders: list[int], levels: list[int]) -> dict[str, object]:
    validate_physical_generators(n)
    parent_caps = root_caps(n)
    all_generators = edge_toggles(n, range(n))
    child_caps = {m: root_caps(m) for m in child_orders}
    projections = {m: projected_masks(n, m) for m in child_orders}
    records: list[dict[str, object]] = []

    for level in levels:
        shell = np.flatnonzero(parent_caps == level)
        target_cap = level - 2
        target_size = int(np.sum(parent_caps <= target_cap))
        arbitrary_distance = quotient_distances(
            parent_caps, target_cap, all_generators
        )
        shell_distance = arbitrary_distance[shell]
        for m in child_orders:
            normalized_level = level / n**1.5
            bad = shell[
                child_caps[m][projections[m][shell]] / m**1.5
                > normalized_level + 1e-12
            ]
            internal_distance = quotient_distances(
                parent_caps,
                target_cap,
                edge_toggles(n, range(m)),
            )
            bad_arbitrary = arbitrary_distance[bad]
            bad_internal = internal_distance[bad]
            records.append(
                {
                    "N": n,
                    "m": m,
                    "cap_level": level,
                    "target_cap": target_cap,
                    "lower_layer_root_size": target_size,
                    "root_shell_size": int(len(shell)),
                    "bad_incidence_count": int(len(bad)),
                    "bad_incidence_mass": float(len(bad) / max(1, len(shell))),
                    "arbitrary_edge": {
                        "distance_histogram_on_shell": distance_histogram(
                            shell_distance
                        ),
                        "distance_histogram_on_bad_incidences": distance_histogram(
                            bad_arbitrary
                        ),
                        "maximum_finite_distance_on_bad_incidences": finite_radius(
                            bad_arbitrary
                        ),
                        "unreachable_bad_incidence_count": int(
                            np.sum(bad_arbitrary == UNREACHABLE)
                        ),
                    },
                    "restricted_internal_edge": {
                        "distance_histogram_on_bad_incidences": distance_histogram(
                            bad_internal
                        ),
                        "maximum_finite_distance_on_bad_incidences": finite_radius(
                            bad_internal
                        ),
                        "unreachable_bad_incidence_count": int(
                            np.sum(bad_internal == UNREACHABLE)
                        ),
                    },
                }
            )

    return {
        "schema": "quadratic-signing-exposed-shell-repair-distance-audit-v1",
        "classification": (
            "exact exhaustive switching-quotient enumeration; distances count "
            "physical edge flips modulo switching"
        ),
        "root_gauge_state_count": int(len(parent_caps)),
        "physical_edge_generator_count": int(len(all_generators)),
        "records": records,
    }


def validate_default_order_eight(payload: dict[str, object]) -> None:
    """Regression-check the exact shell-distance and bottom-layer facts."""
    assert payload["root_gauge_state_count"] == 1 << 21
    assert payload["physical_edge_generator_count"] == math.comb(8, 2)
    records = payload["records"]
    assert isinstance(records, list)
    by_key = {(row["cap_level"], row["m"]): row for row in records}
    shell_histogram = {"1": 97_440, "2": 168_840, "3": 70_560}
    for m in (3, 4, 5):
        bottom = by_key[(10, m)]
        assert bottom["lower_layer_root_size"] == 0
        assert bottom["arbitrary_edge"]["distance_histogram_on_shell"] == {
            "unreachable": 4_200
        }
        exposed = by_key[(12, m)]
        assert (
            exposed["arbitrary_edge"]["distance_histogram_on_shell"]
            == shell_histogram
        )
        assert (
            exposed["arbitrary_edge"]["maximum_finite_distance_on_bad_incidences"]
            == 3
        )
        assert exposed["arbitrary_edge"]["unreachable_bad_incidence_count"] == 0


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--order", type=int, default=8)
    parser.add_argument("--children", type=int, nargs="+", default=[3, 4, 5])
    parser.add_argument("--levels", type=int, nargs="+", default=[10, 12])
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = audit(args.order, args.children, args.levels)
    if (
        args.order == 8
        and args.children == [3, 4, 5]
        and args.levels == [10, 12]
    ):
        validate_default_order_eight(payload)
    rendered = json.dumps(payload, indent=2, sort_keys=True)
    if args.output:
        args.output.write_text(rendered + "\n", encoding="utf-8")
    else:
        print(rendered)


if __name__ == "__main__":
    main()
