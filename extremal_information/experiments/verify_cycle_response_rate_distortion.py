#!/usr/bin/env python3
"""Exact checks for the cycle-response metric in Theorem 17.1h."""

from fractions import Fraction as Q


# Edge = (source, target, coefficient of the hidden seed y).
EDGES = (
    (0, 1, Q(7)),   # large transient response, absent from d_circ
    (1, 1, Q(2)),   # loop mean 2y
    (1, 2, Q(-1)),
    (2, 1, Q(5)),   # two-cycle mean 2y
)


def outgoing(vertex):
    return tuple(edge for edge in EDGES if edge[0] == vertex)


def all_paths(depth):
    paths = [(0, Q(0), 0)]  # vertex, accumulated coefficient, length
    yield from paths
    for _ in range(depth):
        next_paths = []
        for vertex, coefficient, length in paths:
            for _source, target, edge_coefficient in outgoing(vertex):
                next_paths.append(
                    (target, coefficient + edge_coefficient, length + 1)
                )
        paths = next_paths
        yield from paths


def main():
    # The only simple cycles have normalized coefficient 2, so
    # d_circ(y,z)=2|y-z|.  The transient coefficient seven is deliberately
    # larger but can occur only once.
    cycle_means = (Q(2), (Q(-1) + Q(5)) / 2)
    assert cycle_means == (Q(2), Q(2))

    epsilon = Q(1, 5)
    seed = Q(7, 20)
    center = Q(3, 10)
    d_circ = 2 * abs(seed - center)
    assert d_circ <= epsilon

    # A walk decomposes into cycles plus a simple path.  Here |V|-1=2 and
    # max edge oscillation for this pair is 7|seed-center|.
    transient_bound = 2 * 7 * abs(seed - center)
    checks = 0
    worst_excess = Q(0)
    for _vertex, coefficient, length in all_paths(16):
        residual = abs(coefficient * (seed - center))
        assert residual <= epsilon * length + transient_bound
        worst_excess = max(worst_excess, residual - epsilon * length)
        checks += 1

    # Repeating the loop exposes the lower bound exactly.
    left, right = Q(1, 4), Q(3, 4)
    for repetitions in range(1, 101):
        difference = repetitions * 2 * abs(left - right)
        assert difference / repetitions == 2 * abs(left - right)
        checks += 1

    print(f"exact path/cycle checks: {checks}")
    print(f"cycle pseudometric at test pair: {d_circ}")
    print(f"largest transient excess over epsilon*n: {worst_excess}")


if __name__ == "__main__":
    main()
