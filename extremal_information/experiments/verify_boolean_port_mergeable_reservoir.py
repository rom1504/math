#!/usr/bin/env python3
"""Exact finite checks for the min-priority Boolean-port reservoir."""

from collections import Counter
from itertools import permutations, product


def winner(order, subset):
    rank = {u: i for i, u in enumerate(order)}
    return min(subset, key=rank.__getitem__)


def merge_entry(order, left, right):
    rank = {u: i for i, u in enumerate(order)}
    return min((left, right), key=rank.__getitem__)


def main():
    universe = tuple(range(4))
    orders = list(permutations(universe))
    subsets = [
        frozenset(u for u in universe if mask & (1 << u))
        for mask in range(1, 1 << len(universe))
    ]

    # A minimum under a uniform order is exactly uniform on every subset.
    for subset in subsets:
        counts = Counter(winner(order, subset) for order in orders)
        assert set(counts) == set(subset)
        assert len(set(counts.values())) == 1

    # Independent orders give the product law for two replicas.
    for subset in subsets:
        counts = Counter(
            (winner(first, subset), winner(second, subset))
            for first, second in product(orders, repeat=2)
        )
        expected = len(orders) ** 2 // len(subset) ** 2
        assert set(counts) == set(product(subset, repeat=2))
        assert set(counts.values()) == {expected}

    # Coordinatewise minima realize union and are associative.
    for order in orders:
        for left in subsets:
            for right in subsets:
                wl = winner(order, left)
                wr = winner(order, right)
                assert merge_entry(order, wl, wr) == winner(order, left | right)
        for a, b, c in product(universe, repeat=3):
            lhs = merge_entry(order, merge_entry(order, a, b), c)
            rhs = merge_entry(order, a, merge_entry(order, b, c))
            assert lhs == rhs

    print(
        "PASS: exact winner laws and merge associativity on "
        f"{len(subsets)} subsets and {len(orders)} orders"
    )


if __name__ == "__main__":
    main()
