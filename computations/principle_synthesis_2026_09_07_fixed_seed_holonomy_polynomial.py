#!/usr/bin/env python3
"""Exact small-macro raw-kernel contractions from actual Hadamard row laws.

Enumerates 4-regular subgraphs, not full parent spins or all incidence words.
All contraction coefficients and sign comparisons are exact integers/Fractions.
"""
import itertools
import json
import math
from fractions import Fraction

import networkx as nx
import numpy as np


def regular_graphs(n, degree):
    remaining = [degree] * n
    chosen = []

    def visit(i):
        if i == n:
            yield tuple(chosen)
            return
        need = remaining[i]
        available = [j for j in range(i + 1, n) if remaining[j] > 0]
        if need < 0 or need > len(available):
            return
        for neighbors in itertools.combinations(available, need):
            for j in neighbors:
                remaining[j] -= 1
                chosen.append((i, j))
            feasible = all(0 <= remaining[j] <= n - i - 2
                           for j in range(i + 1, n))
            if feasible:
                yield from visit(i + 1)
            for j in reversed(neighbors):
                chosen.pop()
                remaining[j] += 1

    yield from visit(0)


def signing_from_atlas(index):
    graph = nx.graph_atlas_g()[index]
    assert graph.number_of_nodes() == 7
    a = np.ones((8, 8), dtype=np.int64)
    np.fill_diagonal(a, 0)
    for i, j in graph.edges():
        a[i + 1, j + 1] = a[j + 1, i + 1] = -1
    return a


def cap(a):
    spins = np.asarray([(1,) + s for s in
                        itertools.product((-1, 1), repeat=7)], dtype=np.int64)
    return int(np.max(np.abs(np.sum(spins * (spins @ a), axis=1)))) // 2


def row_law():
    h = np.ones((1, 1), dtype=np.int64)
    for _ in range(3):
        h = np.block([[h, h], [h, -h]])
    counts = {}
    total = 0
    for support in itertools.combinations(range(8), 4):
        for signs in itertools.product((-1, 1), repeat=4):
            total += 1
            f = np.zeros(8, dtype=np.int64)
            f[list(support)] = signs
            v = h.T @ f
            if np.all(np.abs(v) == 2):
                key = int(np.sum(v > 0))
                counts[key] = counts.get(key, 0) + 1
    flat = sum(counts.values())
    moments = {}
    for d in range(8):
        value = Fraction()
        for positive, count in counts.items():
            coeff = sum((-1)**j * math.comb(8-positive, j)
                        * math.comb(positive, d-j)
                        for j in range(max(0, d-positive), min(d, 8-positive)+1))
            value += Fraction(count * coeff, flat * math.comb(8, d))
        moments[d] = value
    assert moments == {0: Fraction(1), 1: Fraction(0), 2: Fraction(0),
                       3: Fraction(0), 4: Fraction(-1,35), 5: Fraction(0),
                       6: Fraction(0), 7: Fraction(0)}
    return dict(total_inputs=total, flat_inputs=flat,
                flat_probability=str(Fraction(flat,total)),
                flat_positive_counts=counts,
                flat_moments={d: str(v) for d,v in moments.items()})


def bernstein_coefficients(power_coeffs, left, right):
    """Exact Bernstein coefficients of a polynomial on [left,right]."""
    n = len(power_coeffs) - 1
    translated = [sum(Fraction(power_coeffs[j]) * math.comb(j,k)
                      * left**(j-k) * (right-left)**k
                      for j in range(k,n+1)) for k in range(n+1)]
    return [sum(translated[k] * Fraction(math.comb(i,k),math.comb(n,k))
                for k in range(i+1)) for i in range(n+1)]


def main():
    labels = ["all_positive", "atlas580", "atlas722", "atlas870", "atlas1005"]
    matrices = [np.ones((8,8),dtype=np.int64)-np.eye(8,dtype=np.int64)]
    matrices += [signing_from_atlas(i) for i in [580,722,870,1005]]
    coeffs = {label: {} for label in labels}
    totals = {}
    for s in range(5,9):
        regular = list(regular_graphs(s,4))
        total = 0
        for vertices in itertools.combinations(range(8),s):
            for graph in regular:
                total += 1
                for label, a in zip(labels, matrices):
                    value = math.prod(int(a[vertices[i],vertices[j]]) for i,j in graph)
                    coeffs[label][s] = coeffs[label].get(s,0) + value
        totals[s] = dict(regular_graphs_on_labeled_s=len(regular),
                         embedded_graphs=total)
    results = []
    for label, a in zip(labels,matrices):
        # R(z)-1=z^5 [A5+A6 z+A7 z²+A8 z³].
        polynomial = [coeffs[label][s] for s in range(5,9)]
        bernstein = bernstein_coefficients(polynomial,Fraction(-1,35),Fraction(0))
        if min(bernstein)>0:
            flat_sign = "strictly_below_one_for_all_positive_t"
        elif max(bernstein)<0:
            flat_sign = "strictly_above_one_for_all_positive_t"
        else:
            flat_sign = "not_decided_by_single_Bernstein_interval"
        results.append(dict(label=label, cap=cap(a), coefficients=coeffs[label],
                            flat_interval_bernstein=[str(x) for x in bernstein],
                            flat_contraction_ratio_sign=flat_sign,
                            matrix=a.tolist()))
    output=dict(status="PASS", scope="exact finite row-sector contraction; no asymptotic transfer",
                row_law=row_law(), graph_counts=totals, results=results)
    print(json.dumps(output,indent=2,sort_keys=True))


if __name__ == "__main__":
    main()
