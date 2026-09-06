"""Exact finite checks for the characteristic-three incidence theorem."""

import itertools

import numpy as np


def projective_points(dimension):
    result = []
    for point in itertools.product(range(3), repeat=dimension):
        if any(point) and next(x for x in point if x) == 1:
            result.append(point)
    return np.array(result, dtype=np.int64)


def profile_weights():
    points = np.array(list(itertools.product(range(3), repeat=2)))
    reps = projective_points(2)
    dots = points.dot(reps.T) % 3
    weights = set()
    for mask in range(256):
        spin = np.array([1] + [1 - 2 * ((mask >> i) & 1) for i in range(8)])
        row = [int(spin.sum()) ** 2]
        for j in range(4):
            z = [int(spin[dots[:, j] == value].sum()) for value in range(3)]
            row.append(2 * (sum(x * x for x in z)
                            - z[0] * z[1] - z[0] * z[2] - z[1] * z[2]))
        assert sum(row) == 81
        weights.add(tuple(row))
    weights = np.array(sorted(weights), dtype=np.int64)
    assert len(weights) == 24
    for signs in itertools.product((-1, 1), repeat=4):
        spectrum = np.array((0,) + signs)
        norm_numerator = int(np.max(np.abs(weights.dot(spectrum))))
        assert norm_numerator == (80 if len(set(signs)) == 1 else 72)
    balanced = weights[weights[:, 0] == 1]
    assert len(balanced) == 10
    assert np.array_equal(balanced[:, 1:].sum(axis=0), np.full(4, 200))
    print("PASS F3^2: 24 exact weight profiles; mixed spectra 72/81, monochromatic 80/81")


def plane_incidence():
    points = projective_points(3)
    # In a projective plane every line is the zero set of one dual point.
    incidence = points.dot(points.T) % 3 == 0
    assert incidence.shape == (13, 13)
    assert np.all(incidence.sum(axis=0) == 4)
    assert np.all(incidence.sum(axis=1) == 4)
    for i in range(13):
        for j in range(i):
            assert int(np.sum(incidence[:, i] & incidence[:, j])) == 1
    min_non_odd = 13
    split_pairs = set()
    for mask in range(8192):
        negatives = np.array([(mask >> i) & 1 for i in range(13)], dtype=np.int64)
        counts = incidence.dot(negatives)
        odd_lines = int(np.sum(counts % 2))
        assert odd_lines % 2 == 0
        assert odd_lines <= 12
        min_non_odd = min(min_non_odd, 13 - odd_lines)
        mono = int(np.sum((counts == 0) | (counts == 4)))
        balanced = int(np.sum(counts == 2))
        assert 3 * mono + balanced >= 3
        if mono == 0:
            assert balanced == 3
            assert int(negatives.sum()) in (6, 7)
        split_pairs.add((mono, balanced))
    assert min_non_odd == 1
    assert (0, 3) in split_pairs and (1, 0) in split_pairs
    print("PASS PG(2,3): all 8192 colorings obey parity and sharp 3*N0+N2 >= 3")


def local_law():
    samples = []
    for minority in range(4):
        point = np.ones(4, dtype=np.int64)
        point[minority] = -1
        samples.extend([point, -point])
    samples = np.array(samples)
    assert np.all(samples.sum(axis=0) == 0)
    assert np.array_equal(samples.T.dot(samples), 8 * np.eye(4, dtype=np.int64))
    print("PASS local 3/1 law: exact zero means and identity covariance")


if __name__ == "__main__":
    profile_weights()
    plane_incidence()
    local_law()
