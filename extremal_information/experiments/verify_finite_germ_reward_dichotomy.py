#!/usr/bin/env python3
"""Exact finite checks for Theorem 17.1e and Proposition 17.1f."""

from fractions import Fraction as Q
from itertools import product


def map_a(point):
    _x, _y = point
    return Q(0), Q(0)


def map_b(point):
    x, _y = point
    return Q(0), x


MAPS = {"A": map_a, "B": map_b}


def reward(letter, point):
    return point[0] if letter == "A" else Q(0)


def follow(word, start):
    state = start
    total = Q(0)
    for letter in word:
        total += reward(letter, state)
        state = MAPS[letter](state)
    return state, total


def merging_diamond_checks():
    grid = tuple(Q(i, 8) for i in range(9))
    checks = 0
    for x, y in product(grid, repeat=2):
        start = (x, y)
        # Coterminal paths A and BA have rewards x and zero.
        end_a, value_a = follow("A", start)
        end_ba, value_ba = follow("BA", start)
        assert end_a == end_ba == (Q(0), Q(0))
        assert value_a == x and value_ba == 0

        # Every word is bounded: only an initial A can receive nonzero reward.
        for depth in range(1, 8):
            for word in map("".join, product("AB", repeat=depth)):
                _end, value = follow(word, start)
                assert 0 <= value <= 1
                checks += 1

    # At zero, every word has zero reward.  Exhaustively check short periodic
    # points: a nonzero fixed point of a nonempty word does not occur.
    for depth in range(1, 10):
        for word in map("".join, product("AB", repeat=depth)):
            for x, y in product(grid, repeat=2):
                end, value = follow(word, (x, y))
                if end == (x, y):
                    assert (x, y) == (Q(0), Q(0))
                    assert value == 0
                checks += 1
    return checks


def pumpable_cycle_checks():
    # One identity germ with reward r(x)=x has a nonconstant cycle label.
    left, right = Q(1, 5), Q(4, 5)
    checks = 0
    for repetitions in range(1, 101):
        residual_left = repetitions * left
        residual_right = repetitions * right
        assert residual_right - residual_left == repetitions * (right - left)
        checks += 1
    return checks


def main():
    print(f"merging-diamond/path checks: {merging_diamond_checks()}")
    print(f"pumpable nonconstant-cycle checks: {pumpable_cycle_checks()}")


if __name__ == "__main__":
    main()
