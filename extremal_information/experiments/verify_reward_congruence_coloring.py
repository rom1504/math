#!/usr/bin/env python3
"""Exhaustive small-graph checks for the reward-congruence coloring reduction."""

from itertools import combinations


def partitions(size):
    labels = [0] * size

    def visit(index, maximum):
        if index == size:
            yield tuple(labels)
            return
        for label in range(maximum + 2):
            labels[index] = label
            yield from visit(index + 1, max(maximum, label))

    if size == 0:
        yield ()
    else:
        labels[0] = 0
        yield from visit(1, 0)


def blocks(labels):
    return tuple(
        tuple(vertex for vertex, value in enumerate(labels) if value == label)
        for label in range(max(labels) + 1)
    )


def independent_partition(labels, edges):
    return all(labels[left] != labels[right] for left, right in edges)


def reward_feasible(labels, edges):
    # Twice D(pi) is the maximum within-block range of one edge coordinate.
    for block in blocks(labels):
        for left, right in edges:
            values = tuple(
                1 if vertex == left else -1 if vertex == right else 0
                for vertex in block
            )
            if max(values) - min(values) > 1:
                return False
    return True


def main():
    graphs = 0
    partition_checks = 0
    for size in range(1, 6):
        possible_edges = tuple(combinations(range(size), 2))
        all_partitions = tuple(partitions(size))
        for mask in range(1 << len(possible_edges)):
            edges = tuple(
                edge for index, edge in enumerate(possible_edges)
                if (mask >> index) & 1
            )
            chromatic = size
            reward_optimum = size
            for labels in all_partitions:
                count = max(labels) + 1
                independent = independent_partition(labels, edges)
                feasible = reward_feasible(labels, edges)
                assert feasible == independent
                if independent:
                    chromatic = min(chromatic, count)
                if feasible:
                    reward_optimum = min(reward_optimum, count)
                partition_checks += 1
            assert reward_optimum == chromatic
            graphs += 1
    print(f"graphs checked through five vertices: {graphs}")
    print(f"partition/response checks: {partition_checks}")
    print("minimum half-error reward quotient equals chromatic number")


if __name__ == "__main__":
    main()
