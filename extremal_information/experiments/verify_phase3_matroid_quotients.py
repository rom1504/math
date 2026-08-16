#!/usr/bin/env python3
"""Exhaustive checks for the phase-3 matroid quotient draft.

The default run does four independent finite checks on the subspace
join-semilattices ``L(F_2^w)``:

1. every fixed-subspace projection is a join homomorphism and every one of
   its fibres has the asserted exact rank width;
2. every triggered-contraction relation is a congruence with the asserted
   classes, count, and exact rank width;
3. all join congruences are enumerated through width three (there are 3,616
   at width three), and every one satisfies the kernel decomposition and
   exact oscillation identity of PMQ.2;
4. every enumerated ambient congruence is recovered by pulling its induced
   zero-separating quotient back along the canonical projection.

Vectors are integers with xor as addition.  A subspace is represented by a
``frozenset`` of its vectors.  The calculation is deterministic and uses no
external packages.
"""

from __future__ import annotations

import argparse
import json
from collections import deque
from itertools import combinations
from pathlib import Path


Space = frozenset[int]
Partition = tuple[int, ...]


def span(vectors: object) -> Space:
    values = {0}
    for vector in vectors:
        vector = int(vector)
        values |= {x ^ vector for x in tuple(values)}
    return frozenset(values)


def all_subspaces(w: int) -> list[Space]:
    spaces: set[Space] = {frozenset({0})}
    nonzero = range(1, 1 << w)
    for size in range(1, w + 1):
        for candidate in combinations(nonzero, size):
            candidate_span = span(candidate)
            if len(candidate_span) == 1 << size:
                spaces.add(candidate_span)
    return sorted(spaces, key=lambda item: (len(item), tuple(sorted(item))))


def dimension(space: Space) -> int:
    return len(space).bit_length() - 1


def is_subspace_of(left: Space, right: Space) -> bool:
    return left <= right


def canonical_partition(groups: list[list[int]], n: int) -> Partition:
    owner = [-1] * n
    next_label = 0
    for group in groups:
        if not group:
            continue
        for index in group:
            owner[index] = next_label
        next_label += 1
    assert all(label >= 0 for label in owner)
    remap: dict[int, int] = {}
    result: list[int] = []
    for label in owner:
        if label not in remap:
            remap[label] = len(remap)
        result.append(remap[label])
    return tuple(result)


def partition_groups(partition: Partition) -> list[list[int]]:
    groups: dict[int, list[int]] = {}
    for index, label in enumerate(partition):
        groups.setdefault(label, []).append(index)
    return list(groups.values())


def congruence_closure(
    initial_pairs: list[tuple[int, int]],
    joins: list[list[int]],
) -> Partition:
    n = len(joins)
    parent = list(range(n))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> bool:
        left_root = find(left)
        right_root = find(right)
        if left_root == right_root:
            return False
        if left_root > right_root:
            left_root, right_root = right_root, left_root
        parent[right_root] = left_root
        return True

    for left, right in initial_pairs:
        union(left, right)

    changed = True
    while changed:
        changed = False
        groups: dict[int, list[int]] = {}
        for index in range(n):
            groups.setdefault(find(index), []).append(index)
        for group in groups.values():
            representative = group[0]
            for other in group[1:]:
                for context in range(n):
                    changed |= union(
                        joins[representative][context],
                        joins[other][context],
                    )

    labels: dict[int, int] = {}
    result: list[int] = []
    for index in range(n):
        root = find(index)
        if root not in labels:
            labels[root] = len(labels)
        result.append(labels[root])
    return tuple(result)


def enumerate_congruences(joins: list[list[int]]) -> set[Partition]:
    """Enumerate every congruence by adjoining one generator pair at a time."""

    n = len(joins)
    identity = tuple(range(n))
    seen = {identity}
    queue: deque[Partition] = deque([identity])
    while queue:
        partition = queue.popleft()
        groups = partition_groups(partition)
        representatives = [group[0] for group in groups]
        old_pairs = [
            (group[0], other)
            for group in groups
            for other in group[1:]
        ]
        for left_index in range(len(representatives)):
            for right_index in range(left_index + 1, len(representatives)):
                candidate = congruence_closure(
                    old_pairs
                    + [
                        (
                            representatives[left_index],
                            representatives[right_index],
                        )
                    ],
                    joins,
                )
                if candidate not in seen:
                    seen.add(candidate)
                    queue.append(candidate)
    return seen


def make_lattice(w: int) -> tuple[list[Space], list[list[int]]]:
    spaces = all_subspaces(w)
    index = {space: position for position, space in enumerate(spaces)}
    joins = [
        [index[span(left | right)] for right in spaces]
        for left in spaces
    ]
    return spaces, joins


def check_fixed_projection(max_w: int) -> list[dict[str, int]]:
    records: list[dict[str, int]] = []
    for w in range(max_w + 1):
        spaces, joins = make_lattice(w)
        index = {space: position for position, space in enumerate(spaces)}
        homomorphism_checks = 0
        fibre_checks = 0
        for W in spaces:
            w_index = index[W]
            d = dimension(W)
            images = {joins[x][w_index] for x in range(len(spaces))}
            expected_image_count = len(all_subspaces(w - d))
            assert len(images) == expected_image_count
            for left in range(len(spaces)):
                for right in range(len(spaces)):
                    projected_join = joins[joins[left][right]][w_index]
                    joined_projections = joins[
                        joins[left][w_index]
                    ][joins[right][w_index]]
                    assert projected_join == joined_projections
                    homomorphism_checks += 1
            for image in images:
                fibre = [
                    x for x in range(len(spaces)) if joins[x][w_index] == image
                ]
                ranks = [dimension(spaces[x]) for x in fibre]
                assert max(ranks) - min(ranks) == d
                fibre_checks += 1
        records.append(
            {
                "w": w,
                "flat_count": len(spaces),
                "fixed_subspace_count": len(spaces),
                "homomorphism_checks": homomorphism_checks,
                "exact_fibre_width_checks": fibre_checks,
            }
        )
    return records


def check_triggered_contractions(max_w: int) -> list[dict[str, int]]:
    records: list[dict[str, int]] = []
    for w in range(1, max_w + 1):
        spaces, joins = make_lattice(w)
        strict_pairs = 0
        compatibility_checks = 0
        fibre_checks = 0
        for x_index, X in enumerate(spaces):
            for y_index, Y in enumerate(spaces):
                if X == Y or not is_subspace_of(X, Y):
                    continue
                strict_pairs += 1
                d = dimension(Y) - dimension(X)
                keys: list[tuple[str, int]] = []
                for u_index, U in enumerate(spaces):
                    if is_subspace_of(X, U):
                        keys.append(("inside", joins[u_index][y_index]))
                    else:
                        keys.append(("singleton", u_index))
                labels: dict[tuple[str, int], int] = {}
                partition_list: list[int] = []
                for key in keys:
                    if key not in labels:
                        labels[key] = len(labels)
                    partition_list.append(labels[key])
                partition = tuple(partition_list)

                for left in range(len(spaces)):
                    for right in range(len(spaces)):
                        if partition[left] != partition[right]:
                            continue
                        for context in range(len(spaces)):
                            assert (
                                partition[joins[left][context]]
                                == partition[joins[right][context]]
                            )
                            compatibility_checks += 1

                above_x = sum(is_subspace_of(X, U) for U in spaces)
                above_y = sum(is_subspace_of(Y, U) for U in spaces)
                assert len(set(partition)) == len(spaces) - above_x + above_y

                maximum_width = 0
                for group in partition_groups(partition):
                    ranks = [dimension(spaces[index]) for index in group]
                    width = max(ranks) - min(ranks)
                    maximum_width = max(maximum_width, width)
                    if len(group) > 1:
                        assert width == d
                        fibre_checks += 1
                assert maximum_width == d
                if dimension(X) > 0:
                    assert partition[0] != partition[x_index]
                if Y != spaces[-1]:
                    assert partition[x_index] == partition[y_index]
                    assert any(
                        partition[joins[x_index][z]] != partition[x_index]
                        for z in range(len(spaces))
                    )
        records.append(
            {
                "w": w,
                "flat_count": len(spaces),
                "strict_trigger_pair_count": strict_pairs,
                "congruence_compatibility_checks": compatibility_checks,
                "nontrivial_exact_fibre_width_checks": fibre_checks,
            }
        )
    return records


def restrict_partition(partition: Partition, indices: list[int]) -> Partition:
    remap: dict[int, int] = {}
    result: list[int] = []
    for index in indices:
        label = partition[index]
        if label not in remap:
            remap[label] = len(remap)
        result.append(remap[label])
    return tuple(result)


def check_kernel_classification(max_w: int) -> list[dict[str, int]]:
    records: list[dict[str, int]] = []
    for w in range(max_w + 1):
        spaces, joins = make_lattice(w)
        congruences = enumerate_congruences(joins)
        index = {space: position for position, space in enumerate(spaces)}
        classification_checks = 0
        oscillation_checks = 0
        converse_checks = 0
        for partition in congruences:
            zero_label = partition[0]
            zero_members = [
                spaces[i] for i, label in enumerate(partition) if label == zero_label
            ]
            W = span(vector for member in zero_members for vector in member)
            w_index = index[W]
            d = dimension(W)
            expected_zero_members = {
                i for i, U in enumerate(spaces) if is_subspace_of(U, W)
            }
            actual_zero_members = {
                i for i, label in enumerate(partition) if label == zero_label
            }
            assert actual_zero_members == expected_zero_members

            # Every canonical projection fibre is contained in one theta class.
            for left in range(len(spaces)):
                for right in range(len(spaces)):
                    if joins[left][w_index] == joins[right][w_index]:
                        assert partition[left] == partition[right]
                    # Conversely theta equivalence is completely visible above W.
                    assert (partition[left] == partition[right]) == (
                        partition[joins[left][w_index]]
                        == partition[joins[right][w_index]]
                    )
                    classification_checks += 1

            interval = [
                i for i, U in enumerate(spaces) if is_subspace_of(W, U)
            ]
            induced = restrict_partition(partition, interval)
            assert induced[interval.index(w_index)] != induced[-1] or w_index == interval[-1]
            # More directly: W is alone in the zero class on [W,V].
            induced_zero = induced[interval.index(w_index)]
            assert sum(label == induced_zero for label in induced) == 1

            interval_position = {original: pos for pos, original in enumerate(interval)}
            for group in partition_groups(partition):
                original_residuals = [w - dimension(spaces[i]) for i in group]
                original_width = max(original_residuals) - min(original_residuals)
                image_indices = sorted({joins[i][w_index] for i in group})
                image_labels = {
                    induced[interval_position[image]] for image in image_indices
                }
                assert len(image_labels) == 1
                quotient_residuals = [w - dimension(spaces[i]) for i in image_indices]
                quotient_width = max(quotient_residuals) - min(quotient_residuals)
                assert original_width == d + quotient_width
                oscillation_checks += 1

            # Pulling the induced congruence back along X -> X+W recovers theta.
            pulled: list[int] = []
            for source in range(len(spaces)):
                image = joins[source][w_index]
                pulled.append(induced[interval_position[image]])
            assert canonical_partition(partition_groups(tuple(pulled)), len(spaces)) == partition
            converse_checks += 1

        records.append(
            {
                "w": w,
                "flat_count": len(spaces),
                "enumerated_join_congruence_count": len(congruences),
                "pairwise_kernel_classification_checks": classification_checks,
                "exact_oscillation_checks": oscillation_checks,
                "pullback_reconstruction_checks": converse_checks,
            }
        )
    return records


def run(
    projection_max_w: int,
    trigger_max_w: int,
    classification_max_w: int,
) -> dict[str, object]:
    return {
        "status": "passed",
        "fixed_projection": check_fixed_projection(projection_max_w),
        "triggered_contraction": check_triggered_contractions(trigger_max_w),
        "kernel_classification": check_kernel_classification(classification_max_w),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--projection-max-w", type=int, default=4)
    parser.add_argument("--trigger-max-w", type=int, default=3)
    parser.add_argument("--classification-max-w", type=int, default=3)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name(
            "phase3_matroid_quotients_results.json"
        ),
    )
    args = parser.parse_args()
    result = run(
        args.projection_max_w,
        args.trigger_max_w,
        args.classification_max_w,
    )
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
