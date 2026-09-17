"""Exact Walsh-affine counterexample to an overly strong entropy remainder."""

from __future__ import annotations

from fractions import Fraction
import itertools
import json

import numpy as np


def parity(value):
    return bin(int(value)).count("1") % 2


def graph_words(r):
    side = 2**r
    for entries in itertools.product((0, 1), repeat=r*r):
        matrix_rows = [sum(entries[i*r+j] << j for j in range(r)) for i in range(r)]
        word = np.ones(side*side, dtype=np.int64)
        for u in range(side):
            image = sum(parity(row & u) << i for i, row in enumerate(matrix_rows))
            word[u*side+image] = -1
        yield word


def main():
    reports = []
    for r in range(1, 4):
        side = 2**r
        n = side*side
        walsh = np.asarray([[(-1)**parity(i & j) for j in range(n)]
                            for i in range(n)], dtype=np.int64)
        assert np.array_equal(walsh.T @ walsh, n*np.eye(n, dtype=np.int64))
        words = set()
        responses = set()
        for word in graph_words(r):
            words.add(tuple(word))
            response = Fraction(int(np.abs(walsh @ word).sum()), n)
            assert response == Fraction(3)-Fraction(4, side)
            responses.add(response)
        assert len(words) == 2**(r*r)
        reports.append({"r": r, "physical_dimension": n, "distinct_words": len(words),
                        "exact_mean_absolute_response": str(responses.pop()),
                        "log2_cardinality": r*r})
    r, slab_count = 2, 2
    slab_words = list(graph_words(r))
    n = slab_count*4**r
    walsh = np.asarray([[(-1)**parity(i & j) for j in range(n)]
                        for i in range(n)], dtype=np.int64)
    words = set()
    responses = set()
    for parts in itertools.product(slab_words, repeat=slab_count):
        word = np.concatenate(parts)
        words.add(tuple(word))
        response = Fraction(int(np.abs(walsh @ word).sum()), n)
        assert response <= 1+2*slab_count
        responses.add(response)
    assert len(words) == 2**(slab_count*r*r)
    slab_report = {"r": r, "slabs": slab_count, "physical_dimension": n,
                   "distinct_words": len(words),
                   "responses": sorted(str(x) for x in responses),
                   "proven_response_bound": 1+2*slab_count}
    print(json.dumps({"status": "all exact Fourier and counting checks passed",
                      "cases": reports, "slab_amplification": slab_report}, indent=2))


if __name__ == "__main__":
    main()
