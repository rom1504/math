"""Bounded exact-cap-preserving edge-change test of universal BP-exact.

No claim about all minimizing classes or asymptotic orders. Global parent
minimality is imported from the committed certificates, not solved here.
Every reported cap and child failure is exhaustively integer replayed.
"""

import itertools
import json
from functools import lru_cache
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MINIMUM = {4: 4, 5: 4, 6: 5}


@lru_cache(None)
def spins(n):
    return np.array([(1,) + x for x in itertools.product((-1, 1), repeat=n - 1)], dtype=np.int64)


def cap(a):
    x = spins(len(a))
    return int(np.max(np.abs(np.einsum("bi,ij,bj->b", x, a, x)))) // 2


def balanced_count(a):
    n = len(a)
    m = n // 2
    both = 0
    for left in itertools.combinations(range(n), m):
        if 2 * m == n and 0 not in left:
            continue
        right = tuple(i for i in range(n) if i not in left)
        if cap(a[np.ix_(left, left)]) != MINIMUM[m]:
            continue
        if cap(a[np.ix_(right, right)]) == MINIMUM[n - m]:
            both += 1
    return both


for n in range(3, 9):
    catalog = json.loads((ROOT / f"computations/results/m{n}_minimizer_orbits.json").read_text())
    for item in catalog["classes"]:
        a = np.array(item["representative_matrix"], dtype=np.int64)
        x = spins(n)
        values = np.einsum("bi,ij,bj->b", x, a, x) // 2
        p, r = int(np.max(values)), int(-np.min(values))
        assert max(p, r) == catalog["target_cap"]
        print("catalog", n, "class", item["class"], "oriented_extrema", (p, r), "gap", abs(p - r))

for n, path, key, target in (
    (9, "computations/results/exact_m9.json", "matrix", 12),
    (10, "computations/results/exact_m10.json", "matrix", 13),
    (11, "computations/results/nested_10_in_11_cap17.json", "matrix", 17),
    (12, "computations/results/extension_nested_m11_to_12.json", "parent_matrix", 18),
):
    original = np.array(json.loads((ROOT / path).read_text())[key], dtype=np.int64)
    assert cap(original) == target
    edge_pairs = list(itertools.combinations(range(n), 2))
    x = spins(n)
    edge_terms = np.array([original[i, j] * x[:, i] * x[:, j] for i, j in edge_pairs])
    energies = np.sum(edge_terms, axis=0)
    checked = retained = 0
    minimum_count = balanced_count(original)
    histogram = {}
    found = False
    maximum_orientation_gap = int(abs(np.max(energies) + np.min(energies)))
    for changes in (1, 2, 3):
        for chosen in itertools.combinations(range(len(edge_pairs)), changes):
            changed_energies = energies - 2 * np.sum(edge_terms[list(chosen)], axis=0)
            checked += 1
            if int(np.max(np.abs(changed_energies))) != target:
                continue
            retained += 1
            orientation_gap = int(abs(np.max(changed_energies) + np.min(changed_energies)))
            maximum_orientation_gap = max(maximum_orientation_gap, orientation_gap)
            modified = original.copy()
            for index in chosen:
                i, j = edge_pairs[index]
                modified[i, j] *= -1
                modified[j, i] *= -1
            count = balanced_count(modified)
            histogram[count] = histogram.get(count, 0) + 1
            minimum_count = min(minimum_count, count)
            if count == 0 or orientation_gap >= 4:
                assert cap(modified) == target
                print("FINITE COUNTEREXAMPLE", n, "BP_exact_count", count,
                      "oriented_gap", orientation_gap,
                      "changed_edges", [edge_pairs[index] for index in chosen],
                      "matrix", modified.tolist())
                found = True
                break
        if found:
            break
    print("n", n, "candidates", checked, "retained_exact_cap", retained,
          "minimum_balanced_exact_pairs", minimum_count, "histogram", histogram,
          "maximum_oriented_gap", maximum_orientation_gap)
print("Finite local probe only; no asymptotic or existential-parent falsification.")
