#!/usr/bin/env python3
"""Exact finite witness for the radial optimizer-information ceiling.

The two matrices are the two switching/permutation classes in the certified
order-eight minimizer classification.  This script verifies, using integer
arithmetic except for no operations at all requiring floating point, that

* their complete energy histograms agree;
* their zero-temperature Gibbs overlap second moments differ; and
* their responses to every one-vertex sign field differ.

Run from the repository root with ``python3``.
"""

from collections import Counter
from fractions import Fraction
from itertools import product


A0 = (
    (0, 1, 1, 1, 1, 1, 1, 1),
    (1, 0, 1, 1, 1, 1, -1, -1),
    (1, 1, 0, 1, -1, 1, -1, 1),
    (1, 1, 1, 0, 1, -1, -1, 1),
    (1, 1, -1, 1, 0, -1, 1, -1),
    (1, 1, 1, -1, -1, 0, 1, -1),
    (1, -1, -1, -1, 1, 1, 0, 1),
    (1, -1, 1, 1, -1, -1, 1, 0),
)

A1 = (
    (0, 1, 1, 1, 1, 1, 1, 1),
    (1, 0, 1, 1, -1, -1, 1, -1),
    (1, 1, 0, -1, 1, -1, 1, -1),
    (1, 1, -1, 0, 1, -1, -1, 1),
    (1, -1, 1, 1, 0, -1, -1, 1),
    (1, -1, -1, -1, -1, 0, 1, 1),
    (1, 1, 1, -1, -1, 1, 0, 1),
    (1, -1, -1, 1, 1, 1, 1, 0),
)


def spins(n):
    """One representative of each global-spin pair."""
    return [(1,) + tail for tail in product((-1, 1), repeat=n - 1)]


def energy(A, x):
    return sum(A[i][j] * x[i] * x[j]
               for i in range(len(A)) for j in range(i + 1, len(A)))


def gram(states):
    n = len(states[0])
    return [[sum(x[i] * x[j] for x in states) for j in range(n)]
            for i in range(n)]


def trace_covariance_square(states):
    G = gram(states)
    g = len(states)
    return Fraction(sum(v * v for row in G for v in row), g * g)


def field_response_histogram(A):
    X = spins(len(A))
    H = [energy(A, x) for x in X]
    out = Counter()
    for b in product((-1, 1), repeat=len(A)):
        response = max(
            max(abs(h + field), abs(h - field))
            for h, x in zip(H, X)
            for field in [sum(bi * xi for bi, xi in zip(b, x))]
        )
        out[response] += 1
    return out


def summary(A):
    X = spins(len(A))
    H = [energy(A, x) for x in X]
    ground = [x for h, x in zip(H, X) if abs(h) == max(map(abs, H))]
    return {
        "energy_histogram": dict(sorted(Counter(H).items())),
        "absolute_energy_histogram": dict(
            sorted(Counter(map(abs, H)).items())),
        "cap": max(map(abs, H)),
        "ground_count_projective": len(ground),
        "ground_trace_covariance_square": trace_covariance_square(ground),
        "one_vertex_field_response": dict(sorted(field_response_histogram(A).items())),
    }


def main():
    s0, s1 = summary(A0), summary(A1)
    assert s0["energy_histogram"] == s1["energy_histogram"]
    assert s0["cap"] == s1["cap"] == 10
    assert s0["ground_trace_covariance_square"] == 14
    assert s1["ground_trace_covariance_square"] == 10
    assert s0["one_vertex_field_response"] == {12: 24, 14: 112, 16: 104, 18: 16}
    assert s1["one_vertex_field_response"] == {12: 8, 14: 112, 16: 120, 18: 16}
    print("common absolute-energy histogram:", s0["absolute_energy_histogram"])
    print("ground overlap E[<x,x'>^2]:", 14, 10)
    print("A0 one-vertex response:", s0["one_vertex_field_response"])
    print("A1 one-vertex response:", s1["one_vertex_field_response"])


if __name__ == "__main__":
    main()
