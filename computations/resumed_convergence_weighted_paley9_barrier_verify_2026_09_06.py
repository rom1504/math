"""Exact Paley9 block certificate for the weighted Cayley obstruction."""

import itertools

import numpy as np


points = np.array(list(itertools.product(range(3), repeat=2)), dtype=np.int64)
reps = np.array([[0, 1], [1, 0], [1, 1], [1, 2]], dtype=np.int64)
signs = np.array([1, 1, -1, -1], dtype=np.int64)
raw = np.where(points.dot(reps.T) % 3 == 0, 2, -1).dot(signs)
assert np.all(raw % 3 == 0)
kernel = raw // 3
assert kernel[0] == 0 and np.all(np.abs(kernel[1:]) == 1)
index = {tuple(point): i for i, point in enumerate(points)}
matrix = np.array([[kernel[index[tuple((x-y) % 3)]] for y in points] for x in points])
assert np.array_equal(matrix, matrix.T)
assert np.all(np.diag(matrix) == 0)
assert np.all(matrix.sum(axis=1) == 0)
assert np.array_equal(matrix.dot(matrix), 9*np.eye(9, dtype=np.int64)-np.ones((9, 9), dtype=np.int64))
assert int(np.sum(matrix*matrix)) == 72
energies = []
for mask in range(256):
    spin = np.array([1]+[1-2*((mask >> i) & 1) for i in range(8)], dtype=np.int64)
    energies.append(int(spin.dot(matrix).dot(spin)))
assert max(energies) == 24 and min(energies) == -24
print("PASS exact Paley9: hollow signs, A^2=9I-J, Frobenius squared72, both half-caps12")
for copies in (1, 3, 9, 27):
    order = 9*copies
    frobenius_squared = 72*copies
    cap = 12*copies
    # Squared ratio Q^2/(q*||C||_F^2) equals 2/9 exactly.
    assert 9*cap*cap == 2*order*frobenius_squared
print("PASS exact scalable ratio: Q/(sqrt(q)*Frobenius)=sqrt(2)/3 < 1/2")
