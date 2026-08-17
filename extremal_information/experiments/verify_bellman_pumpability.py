#!/usr/bin/env python3
"""Exact scalar checks for Theorem 17.2b's Bellman/pumping cases."""

from fractions import Fraction as Q


def iterate(lam, defect, start, count):
    value = start
    values = [value]
    for _ in range(count):
        value = defect + lam * value
        values.append(value)
    return values


def main():
    checks = 0

    # Positive defect on a contracting cycle is bounded, not pumpable.
    stable = iterate(Q(1, 2), Q(1), Q(0), 40)
    assert all(value < 2 for value in stable[:-1])
    assert stable[-1] == 2 - Q(1, 2) ** 39
    checks += len(stable)

    # A unit-gain positive defect pumps linearly.
    linear = iterate(Q(1), Q(3, 5), Q(2), 40)
    for index, value in enumerate(linear):
        assert value == 2 + index * Q(3, 5)
        checks += 1

    # Expansion pumps a positive inherited radius even with zero defect.
    exponential = iterate(Q(2), Q(0), Q(1), 30)
    for index, value in enumerate(exponential):
        assert value == 2 ** index
        checks += 1

    # Local absolute data cannot see signed cancellation.
    cancelled = Q(0)
    drifting = Q(0)
    for _ in range(100):
        cancelled += 1
        cancelled -= 1
        drifting += 1
        drifting += 1
        assert cancelled == 0
        checks += 1
    assert drifting == 200

    print(f"exact Bellman/pump/cancellation checks: {checks}")
    print(f"contracting fixed-radius limit: 2 (last={stable[-1]})")
    print(f"unit-gain drift after 100 cycles: {drifting}")


if __name__ == "__main__":
    main()
