#!/usr/bin/env python3
"""Exact checks for exponentially exposed selector-germ families.

The second family uses only three input letters: two permutation selectors
generating the symmetric group and one repeatable scalar probe.
"""

from fractions import Fraction as Q
from itertools import combinations, product
from math import comb


def seed_vector(phases, q):
    values = []
    for phase in phases:
        values.extend(Q(int(index == phase)) for index in range(q))
    return tuple(values)


def rotate_block(vector, block, q):
    result = list(vector)
    start = block * q
    for index in range(q):
        result[start + index] = vector[start + (index - 1) % q]
    return tuple(result)


def maxplus_selector(vector, block, q):
    """Finite matrix: intended selector weight 0, all off entries -2."""

    dimension = len(vector)
    result = []
    for output in range(dimension):
        intended = output
        if output // q == block:
            intended = block * q + (output - block * q - 1) % q
        candidates = tuple(
            vector[source] + (0 if source == intended else -2)
            for source in range(dimension)
        )
        assert candidates[intended] > max(
            value for source, value in enumerate(candidates)
            if source != intended
        )
        result.append(candidates[intended])
    return tuple(result)


def probe(vector, block, coordinate, q):
    start = block * q
    mean = sum(vector[start:start + q], Q(0)) / q
    return vector[start + coordinate] - mean


def verify(q, m):
    phases = tuple(product(range(q), repeat=m))
    vectors = {phase: seed_vector(phase, q) for phase in phases}
    assert len(set(vectors.values())) == q ** m
    checks = 0

    # The finite max-plus matrices exactly realize commuting block rotations.
    for phase, vector in vectors.items():
        for block in range(m):
            expected = rotate_block(vector, block, q)
            assert maxplus_selector(vector, block, q) == expected
            checks += 1
        for left in range(m):
            for right in range(left + 1, m):
                assert rotate_block(
                    rotate_block(vector, left, q), right, q
                ) == rotate_block(
                    rotate_block(vector, right, q), left, q
                )
                checks += 1

    # Every pair is separated by one on a repeatable identity-probe loop.
    for phase, other in combinations(phases, 2):
        vector, other_vector = vectors[phase], vectors[other]
        separating = next(block for block in range(m)
                          if phase[block] != other[block])
        coordinate = phase[separating]
        difference = abs(
            probe(vector, separating, coordinate, q)
            - probe(other_vector, separating, coordinate, q)
        )
        assert difference == 1
        checks += 1

    return len(phases), checks


def main():
    total = 0
    for q in (2, 3):
        for m in range(1, 6):
            germs, checks = verify(q, m)
            assert germs == q ** m
            total += checks
    print(f"exact group/action/probe checks: {total}")
    print("family size: q^m; all pairwise cycle-response distances equal 1")
    print("therefore epsilon<1/2 forces q^m predictor states")

    fixed_checks = 0
    for dimension in range(2, 11):
        weight = dimension // 2
        seeds = tuple(
            tuple(Q(int(index in support)) for index in range(dimension))
            for support in combinations(range(dimension), weight)
        )
        assert len(seeds) == comb(dimension, weight)

        # A long cycle and an adjacent transposition generate S_r.  Breadth
        # first search verifies that their induced orbit is the entire
        # constant-weight layer; the max-plus realization is checked at every
        # reached point.
        cycle = tuple((index - 1) % dimension for index in range(dimension))
        swap = tuple(1 if index == 0 else 0 if index == 1 else index
                     for index in range(dimension))

        def select(vector, permutation):
            result = []
            for output, intended in enumerate(permutation):
                candidates = tuple(
                    vector[source] + (0 if source == intended else -2)
                    for source in range(dimension)
                )
                assert candidates[intended] > max(
                    value for source, value in enumerate(candidates)
                    if source != intended
                )
                result.append(candidates[intended])
            return tuple(result)

        orbit = {seeds[0]}
        frontier = [seeds[0]]
        while frontier:
            vector = frontier.pop()
            for permutation in (cycle, swap):
                image = select(vector, permutation)
                fixed_checks += 1
                if image not in orbit:
                    orbit.add(image)
                    frontier.append(image)
        assert orbit == set(seeds)

        # A coordinate in the symmetric difference can be moved to coordinate
        # zero by a permutation word.  Repeating the identity probe there has
        # reward-rate separation exactly one.
        mean = Q(weight, dimension)
        for left, right in combinations(seeds, 2):
            coordinate = next(
                index for index in range(dimension)
                if left[index] != right[index]
            )
            difference = abs((left[coordinate] - mean)
                             - (right[coordinate] - mean))
            assert difference == 1
            fixed_checks += 1

    print(f"fixed-alphabet orbit/probe checks: {fixed_checks}")
    print("two permutation inputs plus one probe expose binom(r,floor(r/2)) states")


if __name__ == "__main__":
    main()
