"""Exact finite diagnostic; NOT an asymptotic reverse-Fekete falsifier.

Input is the committed original order-14 signing named below. No output
files are written. Fixing spin zero removes only the global reversal.
"""

import itertools
import json
from collections import Counter
from pathlib import Path

import numpy as np


def cap(matrix):
    n = len(matrix)
    tails = np.array(list(itertools.product((-1, 1), repeat=n - 1)), dtype=np.int64)
    spins = np.column_stack((np.ones(len(tails), dtype=np.int64), tails))
    doubled = np.einsum("bi,ij,bj->b", spins, matrix, spins)
    assert np.all(doubled % 2 == 0)
    return int(np.max(np.abs(doubled))) // 2


path = Path(__file__).resolve().parent / "results/heuristic_m14_from_conference.json"
matrix = np.array(json.loads(path.read_text())["matrix"], dtype=np.int64)
assert matrix.shape == (14, 14)
assert np.array_equal(matrix, matrix.T)
assert not np.any(np.diag(matrix))
assert np.all(np.abs(matrix[np.triu_indices(14, 1)]) == 1)
assert cap(matrix) == 21
histogram = Counter()
for tail in itertools.combinations(range(1, 14), 6):
    left = (0,) + tail
    right = tuple(i for i in range(14) if i not in left)
    histogram[cap(matrix[np.ix_(left, left)]), cap(matrix[np.ix_(right, right)])] += 1
assert histogram == {(9, 9): 624, (11, 11): 1092}
left = (0, 8, 9, 10, 11, 12, 13)
right = (1, 2, 3, 4, 5, 6, 7)
assert cap(matrix[np.ix_(left, left)]) == cap(matrix[np.ix_(right, right)]) == 9
# Exact squared/cubed comparisons: child cap below half, powered gap positive.
assert 4 * 9**2 < 7**3
assert 21**2 < 8 * 9**2
print("PASS: parent cap 21; 624 sub-half child pairs (9,9); 1092 pairs (11,11).")
print("Finite-order witness only; no asymptotic o(N)-defect conclusion.")
