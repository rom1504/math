"""Exact ground-code common-mode tests on stored finite minimizers.

Global optimality is imported from the old certificates; caps, all ground
states, fourth moments, and edge-character equivalences are checked anew.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
import importlib.util
import itertools
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def maximum_matching_size(edges):
    reachable = {0: 0}
    for i, j in edges:
        bits = (1 << i) | (1 << j)
        for mask, size in list(reachable.items()):
            if not mask & bits:
                reachable[mask | bits] = max(reachable.get(mask | bits, 0), size+1)
    return max(reachable.values())


def main():
    source = ROOT / "computations/transfer_adversary_minimizer_isotropy_2026_09_06.py"
    spec = importlib.util.spec_from_file_location("stored_minimizers", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    report = []
    for name, matrix, provenance in module.cases():
        n = len(matrix)
        if n % 2:
            continue
        a = np.asarray(matrix, dtype=np.int64)
        words = np.asarray([(1,)+row for row in itertools.product((-1, 1), repeat=n-1)],
                           dtype=np.int64)
        energies = np.einsum("bi,ij,bj->b", words, a, words)//2
        cap = int(np.abs(energies).max())
        assert cap == provenance["minimum_cap_imported"]
        ground = words[np.abs(energies) == cap]
        fourth_sum = max(abs(int(np.prod(ground[:, indices], axis=1).sum()))
                         for indices in itertools.combinations(range(n), 4))
        fourth = Fraction(fourth_sum, len(ground))
        classes = defaultdict(list)
        for i, j in itertools.combinations(range(n), 2):
            column = ground[:, i]*ground[:, j]
            column *= column[0]
            classes[tuple(column)].append((i, j))
        largest_matching = max(maximum_matching_size(edges) for edges in classes.values())
        if n >= 6:
            assert len(classes) == n*(n-1)//2
        report.append({"case": name, "n": n, "cap": cap,
                       "projective_ground_count": len(ground),
                       "max_absolute_fourth_moment": str(fourth),
                       "edge_character_classes_up_to_sign": len(classes),
                       "largest_common_character_matching": largest_matching,
                       "common_perfect_matching_mode_possible": largest_matching == n//2})
    print(json.dumps({"status": "all integer/rational checks passed", "cases": report}, indent=2))


if __name__ == "__main__":
    main()
