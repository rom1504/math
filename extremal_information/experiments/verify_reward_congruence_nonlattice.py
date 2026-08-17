#!/usr/bin/env python3
"""Exact checks for the bounded-reward congruence nonlattice example."""

from fractions import Fraction as Q
from itertools import product


STATES = ("I", "A", "B")
LETTERS = ("a", "b")


def transition(state, letter):
    del state
    return "A" if letter == "a" else "B"


def reward(state, letter):
    return Q(int(state == "A" and letter == "b"))


def raw_total(state, word):
    total = Q(0)
    for letter in word:
        total += reward(state, letter)
        state = transition(state, letter)
    return total


def quotient_total(state, word, partition, toll):
    total = Q(0)
    for letter in word:
        block = partition[state]
        total += toll[(block, letter)]
        state = transition(state, letter)
    return total


def check_partition(partition, toll, bound, max_depth=12):
    # Forward congruence.
    for left in STATES:
        for right in STATES:
            if partition[left] != partition[right]:
                continue
            for letter in LETTERS:
                assert partition[transition(left, letter)] == partition[
                    transition(right, letter)
                ]

    maximum = Q(0)
    for depth in range(max_depth + 1):
        for word in product(LETTERS, repeat=depth):
            for state in STATES:
                error = abs(
                    raw_total(state, word)
                    - quotient_total(state, word, partition, toll)
                )
                assert error <= bound
                maximum = max(maximum, error)
    return maximum


def main():
    left_partition = {"I": 0, "A": 0, "B": 1}
    left_toll = {
        (0, "a"): Q(0), (0, "b"): Q(1),
        (1, "a"): Q(0), (1, "b"): Q(0),
    }
    right_partition = {"I": 1, "A": 0, "B": 1}
    right_toll = {
        (0, "a"): Q(0), (0, "b"): Q(1),
        (1, "a"): Q(0), (1, "b"): Q(0),
    }
    assert check_partition(left_partition, left_toll, Q(1)) == 1
    assert check_partition(right_partition, right_toll, Q(0)) == 0

    # Their join is the one-block partition.  The constant words force both
    # letter tolls to zero; then (ba)^n has raw reward n from A.
    for repetitions in range(1, 13):
        word = tuple("ba" * repetitions)
        assert raw_total("A", word) == repetitions

    # Every pair of raw states has only a one-step response discrepancy,
    # because the first input resets all of them to the same successor.
    pair_bound = Q(0)
    for depth in range(13):
        for word in product(LETTERS, repeat=depth):
            for left in STATES:
                for right in STATES:
                    pair_bound = max(
                        pair_bound,
                        abs(raw_total(left, word) - raw_total(right, word)),
                    )
    assert pair_bound == 1
    print("two incomparable two-block congruences have errors 1 and 0")
    print("their one-block join has linear error on (ba)^n")
    print("all pairwise same-word raw response differences stay at most 1")


if __name__ == "__main__":
    main()
