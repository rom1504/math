#!/usr/bin/env python3
"""Exact toy checks for the selector-cycle periodicity theorem."""

from fractions import Fraction as Q


def affine_iterate(value, slope, shift, depth):
    for _ in range(depth):
        value = slope * value + shift
    return value


def compact_periodic_checks():
    checks = 0

    # A coordinate swap with offset induces z -> 1-z on projective R^2/R1.
    # Its selector has period two and every point of [0,1] is projectively
    # periodic with period dividing two.
    for numerator in range(65):
        z = Q(numerator, 64)
        assert 0 <= z <= 1
        assert affine_iterate(z, Q(-1), Q(1), 2) == z
        assert 0 <= affine_iterate(z, Q(-1), Q(1), 1) <= 1
        checks += 1

    # Unequal means on the two identity-selector cycles induce z -> z+1/5.
    # On the compact domain [0,1], at most six consecutive starts can remain
    # legal; there is no projectively periodic point.
    step = Q(1, 5)
    maximum = 0
    for denominator in range(1, 81):
        for numerator in range(denominator + 1):
            z = Q(numerator, denominator)
            length = 0
            while 0 <= z + length * step <= 1:
                length += 1
            maximum = max(maximum, length)
            checks += 1
    assert maximum == 6
    assert all(affine_iterate(Q(n, 20), Q(1), step, 1) != Q(n, 20)
               for n in range(21))

    # Compactness is indispensable: x -> x-1 on (0,infinity) realizes every
    # finite repetition length from a depth-dependent seed, but no one seed
    # remains legal forever.
    for depth in range(1, 101):
        seed = Q(depth + 1)
        assert all(affine_iterate(seed, Q(1), Q(-1), time) > 0
                   for time in range(depth + 1))
        assert affine_iterate(seed, Q(1), Q(-1), depth + 1) == 0
        checks += 1

    return checks


def functional_graph_mean_checks():
    # sigma=(1,1,2): coordinate zero feeds the fixed cycle at one, while two
    # is another fixed cycle.  Formula (3.8) is checked directly.
    sigma = (1, 1, 2)
    shift = (Q(3), Q(2), Q(2))
    vector = (Q(-4), Q(5), Q(7))

    def update(state):
        return tuple(state[sigma[index]] + shift[index]
                     for index in range(3))

    state = vector
    checks = 0
    history = [state]
    for _ in range(20):
        state = update(state)
        history.append(state)
    # Both terminal selector cycles have mean two.  After the one-step tree
    # transient, every further step adds the common projective gauge 2.
    for time in range(1, 20):
        assert tuple(history[time][i] + 2 for i in range(3)) \
            == history[time + 1]
        checks += 1

    unequal_shift = (Q(3), Q(2), Q(-1))

    def update_unequal(state):
        return tuple(state[sigma[index]] + unequal_shift[index]
                     for index in range(3))

    state = vector
    differences = []
    for _ in range(12):
        differences.append(state[1] - state[2])
        state = update_unequal(state)
    assert all(differences[index + 1] - differences[index] == 3
               for index in range(len(differences) - 1))
    checks += len(differences) - 1
    return checks


def main():
    print(f"compact/noncompact cycle checks: {compact_periodic_checks()}")
    print(f"functional-graph mean checks: {functional_graph_mean_checks()}")


if __name__ == "__main__":
    main()
