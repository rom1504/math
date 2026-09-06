#!/usr/bin/env python3
"""Finite checks for the Wave 33 selector-transport certificate."""

from itertools import combinations, product
from math import comb
import random

from forest_mincut_r32_check import formula_cost


def dpr_partial(y, z, overlap):
    if not overlap:
        return 0
    disagree = sum(y[i] != z[i] for i in overlap)
    return min(disagree, len(overlap) - disagree)


def transport_bound(ncoord, edges, z, masks, words):
    # Node 0 is the full root; selector node v+1 uses masks[v], words[v].
    node_masks = [set(range(ncoord))] + list(masks)
    node_words = [tuple(z)] + list(words)
    total = 0
    for u, v in edges:
        # Test trees are indexed so u is the parent of v.
        overlap = node_masks[u] & node_masks[v]
        total += dpr_partial(node_words[u], node_words[v], overlap)
        total += len(node_masks[v] - node_masks[u])
    return total


def random_transport_checks():
    rng = random.Random(3305)
    for _ in range(240):
        ncoord = rng.randint(1, 5)
        q = rng.randint(1, 4)
        # Random recursive rooted tree, with parent index smaller than child.
        edges = [(rng.randrange(v), v) for v in range(1, q + 1)]
        z = tuple(rng.choice((-1, 1)) for _ in range(ncoord))
        masks = []
        words = []
        for _v in range(q):
            mask = {i for i in range(ncoord) if rng.randrange(2)}
            if not mask:
                mask = {rng.randrange(ncoord)}
            masks.append(mask)
            words.append(tuple(rng.choice((-1, 1)) for _ in range(ncoord)))
        exact = formula_cost(ncoord, edges, z, masks, words)
        bound = transport_bound(ncoord, edges, z, masks, words)
        assert exact <= bound, (ncoord, edges, z, masks, words, exact, bound)


def johnson_neighbors(s, n, m, h):
    return sum(1 for t in combinations(range(n), m)
               if len(set(s) - set(t)) == h)


def sphere_count_checks():
    for n in range(4, 10):
        for m in range(1, n):
            s = tuple(range(m))
            for h in range(min(m, n - m) + 1):
                assert johnson_neighbors(s, n, m, h) == comb(m, h) * comb(n - m, h)


def tuple_count_check():
    # Exhaustively verify (J33.3) for a small slice and s=3.
    n, m, sample_size = 6, 3, 3
    selectors = [set(s) for s in combinations(range(n), m)]
    N = len(selectors)

    def dj(a, b):
        return len(a - b)

    # For three labeled points, the MST is the sum of the two smallest of
    # the three pair distances.
    for D in range(4):
        good = 0
        for triple in product(range(N), repeat=sample_size):
            ds = sorted((dj(selectors[triple[0]], selectors[triple[1]]),
                         dj(selectors[triple[0]], selectors[triple[2]]),
                         dj(selectors[triple[1]], selectors[triple[2]])))
            good += ds[0] + ds[1] <= D
        upper_count = N * sample_size ** (sample_size - 2)
        upper_count *= comb(D + sample_size - 1, sample_size - 1) * n ** (2 * D)
        assert good <= upper_count


if __name__ == "__main__":
    random_transport_checks()
    sphere_count_checks()
    tuple_count_check()
    print("PASS johnson_transport_r33")
