"""Exact two-copy overlap decompositions for connected graph observables.

Uses unrestricted parent choices, exact uniform subsets, and rational
normalization. No output files and no probabilistic simulation.
"""

from fractions import Fraction as F
from functools import lru_cache
import itertools
import json
import math

import numpy as np


def falling(n, k):
    return 0 if k > n else math.factorial(n)//math.factorial(n-k)


def matchings(first, second):
    for size in range(min(first, second)+1):
        for left in itertools.combinations(range(first), size):
            for right in itertools.combinations(range(second), size):
                for image in itertools.permutations(left):
                    yield dict(zip(right, image))


def merged(g, h, matching):
    first, edges_g = g
    second, edges_h = h
    next_label = first
    labels = []
    for j in range(second):
        if j in matching:
            labels.append(matching[j])
        else:
            labels.append(next_label)
            next_label += 1
    edges = tuple(sorted(tuple(sorted(edge)) for edge in
                         list(edges_g)+[(labels[u], labels[v]) for u, v in edges_h]))
    return next_label, edges


def evaluate(parent, graph, subset):
    vertices, edges = graph
    return sum(math.prod(int(parent[colors[u], colors[v]]) for u, v in edges)
               for colors in itertools.permutations(subset, vertices))


def connected(graph):
    vertices, edges = graph
    seen = {0}
    while True:
        expanded = seen | {v for u, v in edges if u in seen} | {u for u, v in edges if v in seen}
        if expanded == seen:
            return len(seen) == vertices
        seen = expanded


def main():
    rng = np.random.default_rng(260906)
    graphs = [
        (1, ()),
        (2, ((0, 1), (0, 1))),
        (3, ((0, 1), (1, 2), (0, 2))),
        (4, ((0, 1), (1, 2), (2, 3), (0, 3))),
        (1, ((0, 0),)),
        (2, ((0, 0), (0, 1), (0, 1))),
        (3, ((0, 1), (0, 1), (0, 2), (0, 2))),
    ]
    ambient_checks = covariance_checks = connectivity_checks = 0
    maximum_covariance = F(0)
    for order in (4, 6):
        for _ in range(4):
            parent = rng.choice([-1, 1], (order, order))
            parent = np.triu(parent)+np.triu(parent, 1).T

            @lru_cache(None)
            def total(graph):
                return evaluate(parent, graph, tuple(range(order)))

            for g, h in itertools.product(graphs, repeat=2):
                k, ell = len(g[1]), len(h[1])
                if (k+ell) % 2:
                    continue
                overlaps = []
                for matching in matchings(g[0], h[0]):
                    union = merged(g, h, matching)
                    if matching:
                        assert connected(union)
                        a_union = F(union[0]-1)-F(k+ell, 2)
                        a_sum = F(g[0]+h[0]-2)-F(k+ell, 2)
                        assert a_union == a_sum-len(matching)+1
                        connectivity_checks += 1
                    overlaps.append((len(matching), union, total(union)))
                # Equation (7) before division by ambient falling factorials.
                assert total(g)*total(h) == sum(value for _, _, value in overlaps)
                ambient_checks += 1
                for size in range(2, order+1):
                    subsets = list(itertools.combinations(range(order), size))
                    values_g = [evaluate(parent, g, subset) for subset in subsets]
                    values_h = [evaluate(parent, h, subset) for subset in subsets]
                    norm = size**(2+(k+ell)//2)
                    direct_product = F(sum(x*y for x, y in zip(values_g, values_h)), len(subsets)*norm)
                    overlap_product = sum((F(falling(size, union[0]), falling(order, union[0]))*value/norm
                                           for _, union, value in overlaps if union[0] <= size), F(0))
                    assert direct_product == overlap_product
                    mean_product = F(sum(values_g)*sum(values_h), len(subsets)**2*norm)
                    covariance = direct_product-mean_product
                    maximum_covariance = max(maximum_covariance, abs(covariance))
                    if g == h:
                        assert covariance >= 0
                    # The covariance decomposition separates actual sample
                    # overlaps from the disjoint-versus-product correction.
                    sample_overlap = sum((F(falling(size, union[0]), falling(order, union[0]))*value/norm
                                          for z, union, value in overlaps if z and union[0] <= size), F(0))
                    disjoint = next((union, value) for z, union, value in overlaps if z == 0)
                    disjoint_term = (F(falling(size, disjoint[0][0]), falling(order, disjoint[0][0]))
                                     *disjoint[1]/norm if disjoint[0][0] <= size else F(0))
                    assert covariance == sample_overlap+disjoint_term-mean_product
                    covariance_checks += 1
    print(json.dumps({"seed": 260906, "ambient_overlap_identities": ambient_checks,
                      "exact_subset_covariance_decompositions": covariance_checks,
                      "connected_overlap_exponent_checks": connectivity_checks,
                      "largest_finite_absolute_covariance": str(maximum_covariance),
                      "status": "all exact two-copy decompositions pass"}, indent=2))


if __name__ == "__main__":
    main()
