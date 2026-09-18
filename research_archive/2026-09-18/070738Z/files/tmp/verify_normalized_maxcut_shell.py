#!/usr/bin/env python3
"""Exact finite checks for the unit-load max-plus shell theorem."""

from itertools import product
import random


def hamming(x, y):
    return sum(a != b for a, b in zip(x, y))


def canonical(x):
    nx = tuple(-a for a in x)
    return min(tuple(x), nx)


def projective_words(w):
    return sorted({canonical(x) for x in product((-1, 1), repeat=w)})


def dproj(x, y, weights):
    direct = sum(c for a, b, c in zip(x, y, weights) if a != b)
    complement = sum(c for a, b, c in zip(x, y, weights) if a == b)
    return min(direct, complement)


def shell_response(f, words, weights, s):
    # Constants C_f+2 sum_i ell_i are omitted.
    oriented = list(product((-1, 1), repeat=len(weights)))
    return max(f[canonical(y)] - sum(c for a, b, c in zip(s, y, weights) if a != b)
               for y in oriented)


def check_trial(w, rng):
    words = projective_words(w)
    weights = tuple(rng.randint(0, 4) for _ in range(w))
    # A random max of weighted distance cones is automatically Lipschitz.
    seeds = [(rng.choice(words), rng.randint(-7, 7)) for _ in range(8)]
    f = {x: max(b - dproj(x, a, weights) for a, b in seeds) for x in words}

    for x, y in product(words, repeat=2):
        assert abs(f[x] - f[y]) <= dproj(x, y, weights)
    for s in product((-1, 1), repeat=w):
        assert shell_response(f, words, weights, s) == f[canonical(s)]

    deltas = []
    for i in range(w):
        delta = 0
        for s in product((-1, 1), repeat=w):
            t = list(s)
            t[i] *= -1
            delta = max(delta, abs(f[canonical(s)] - f[canonical(tuple(t))]))
        deltas.append(delta)
    for x, y in product(words, repeat=2):
        assert abs(f[x] - f[y]) <= dproj(x, y, deltas)


def main():
    rng = random.Random(20260816)
    checks = 0
    for w in range(2, 7):
        for _ in range(200):
            check_trial(w, rng)
            checks += 1
    print(f"PASS: {checks} anisotropic projective shell trials (w=2..6)")


if __name__ == "__main__":
    main()
