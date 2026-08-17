#!/usr/bin/env python3
"""Finite checks for atomic type quantization on signed-sum zonotopes."""

from __future__ import annotations

import itertools
import math
import random
from collections import Counter


def dot(x, y):
    return sum(a * b for a, b in zip(x, y))


def support_formula(vectors, theta):
    return sum(abs(dot(theta, vector)) for vector in vectors)


def support_bruteforce(vectors, theta):
    return max(
        dot(theta, tuple(sum(sign * vector[j] for sign, vector in zip(signs, vectors))
                         for j in range(len(theta))))
        for signs in itertools.product((-1, 1), repeat=len(vectors))
    )


def quantize(vector, spacing):
    return tuple(round(value / spacing) * spacing for value in vector)


def histogram(vectors, spacing):
    return Counter(quantize(vector, spacing) for vector in vectors)


def decoded_support(hist, theta):
    return sum(count * abs(dot(theta, centre)) for centre, count in hist.items())


def random_unit_ball_vector(rng, dimension):
    while True:
        vector = tuple(rng.uniform(-1.0, 1.0) for _ in range(dimension))
        if sum(value * value for value in vector) <= 1.0:
            return vector


def main():
    rng = random.Random(20260816)
    dimension = 2
    spacing = 0.2
    atom_radius = math.sqrt(dimension) * spacing / 2
    directions = []
    for k in range(20):
        angle = 2 * math.pi * k / 20
        directions.append((math.cos(angle), math.sin(angle)))

    support_checks = 0
    error_checks = 0
    merge_checks = 0
    for n in range(1, 11):
        for _ in range(120):
            vectors = tuple(random_unit_ball_vector(rng, dimension) for _ in range(n))
            hist = histogram(vectors, spacing)
            cut = rng.randrange(n + 1)
            merged = histogram(vectors[:cut], spacing) + histogram(vectors[cut:], spacing)
            assert merged == hist
            merge_checks += 1

            for theta in directions:
                formula = support_formula(vectors, theta)
                brute = support_bruteforce(vectors, theta)
                decoded = decoded_support(hist, theta)
                assert abs(formula - brute) <= 1e-10
                assert abs(formula - decoded) <= n * atom_radius + 1e-10
                support_checks += 1
                error_checks += 1

    print(f"signed-support identities: {support_checks}")
    print(f"root-scale atom-net error checks: {error_checks}")
    print(f"histogram merge checks: {merge_checks}")


if __name__ == "__main__":
    main()
