"""Integer-only identities for complement exposure in actual sign weaves."""

from itertools import product
import json

import numpy as np


def hadamard(order):
    matrix = np.ones((1, 1), dtype=np.int64)
    while len(matrix) < order:
        matrix = np.block([[matrix, matrix], [matrix, -matrix]])
    assert len(matrix) == order
    return matrix


def components(gram):
    unseen = set(range(len(gram)))
    count = 0
    while unseen:
        count += 1
        stack = [unseen.pop()]
        while stack:
            vertex = stack.pop()
            adjacent = {j for j in unseen if gram[vertex, j] != 0}
            unseen.difference_update(adjacent)
            stack.extend(adjacent)
    return count


def check(bases, seed, deleted):
    m = len(seed)
    kept = m-deleted
    full = np.empty((m*m, m*m), dtype=np.int64)
    for i in range(m):
        for j in range(m):
            full[i*m:(i+1)*m, j*m:(j+1)*m] = seed[i, j]*np.outer(
                bases[i][:, j], bases[j][:, i])
    assert np.array_equal(full, full.T)
    assert np.array_equal(full@full, m*m*np.eye(m*m, dtype=np.int64))
    retain = [i*m+a for i in range(m) for a in range(kept)]
    remove = [i*m+a for i in range(m) for a in range(kept, m)]
    kernel = full[np.ix_(retain, retain)]
    cross = full[np.ix_(retain, remove)]
    complement = full[np.ix_(remove, remove)]
    epsilon = complement[::deleted, ::deleted]
    retained_columns, deleted_columns = [], []
    gram_components = []
    row_ratios = []
    for i in range(m):
        retained = cross[i*kept:(i+1)*kept, ::deleted]
        normalized_deleted = np.column_stack([
            cross[j*kept, i*deleted:(i+1)*deleted]*cross[j*kept, i*deleted]
            for j in range(m)])
        assert np.all(normalized_deleted[0] == 1)
        retained_columns.append(retained)
        deleted_columns.append(normalized_deleted)
        retained_gram = retained.T@retained
        deleted_gram = normalized_deleted.T@normalized_deleted
        assert np.array_equal(retained_gram + epsilon[i, :, None]*deleted_gram*epsilon[i, None, :],
                              m*np.eye(m, dtype=np.int64))
        count = components(deleted_gram)
        assert count <= deleted
        gram_components.append(count)
        if deleted % 2:
            assert np.all(deleted_gram % 2 == 1)
            assert count == 1
            ratio = np.ones(m, dtype=np.int64)
            for j in range(1, m):
                numerator, denominator = -retained_gram[0, j], deleted_gram[0, j]
                assert numerator % denominator == 0
                ratio[j] = numerator//denominator
            assert np.all(np.abs(ratio) == 1)
            assert np.array_equal(ratio, epsilon[i, 0]*epsilon[i])
            row_ratios.append(ratio)
    reconstructed = np.empty_like(kernel)
    for i in range(m):
        for j in range(m):
            reconstructed[i*kept:(i+1)*kept, j*kept:(j+1)*kept] = (
                epsilon[i, j]*np.outer(retained_columns[i][:, j], retained_columns[j][:, i]))
    assert np.array_equal(reconstructed, kernel)
    if deleted % 2:
        q = np.array(row_ratios)
        recovered_epsilon = q[0, :, None]*q
        assert np.array_equal(recovered_epsilon, epsilon[0, 0]*epsilon)
        recovered = np.empty_like(kernel)
        for i in range(m):
            for j in range(m):
                recovered[i*kept:(i+1)*kept, j*kept:(j+1)*kept] = (
                    recovered_epsilon[i, j]*np.outer(retained_columns[i][:, j], retained_columns[j][:, i]))
        assert np.array_equal(recovered, epsilon[0, 0]*kernel)
    return max(gram_components)


def random_bases(m, rng):
    base = hadamard(m)
    return [rng.choice([-1, 1], m)[:, None]*base[rng.permutation(m)][:, rng.permutation(m)]
            *rng.choice([-1, 1], m)[None, :] for _ in range(m)]


def main():
    rng = np.random.default_rng(260906)
    counts = {}
    max_components = {}
    m = 4
    bases = random_bases(m, rng)
    upper = [(i, j) for i in range(m) for j in range(i, m)]
    for signs in product((-1, 1), repeat=len(upper)):
        seed = np.empty((m, m), dtype=np.int64)
        for (i, j), value in zip(upper, signs):
            seed[i, j] = seed[j, i] = value
        for deleted in (1, 2, 3):
            key = f"m{m}_deleted{deleted}"
            observed = check(bases, seed, deleted)
            counts[key] = counts.get(key, 0)+1
            max_components[key] = max(max_components.get(key, 0), observed)
    for m in (8, 16):
        for _ in range(32):
            bases = random_bases(m, rng)
            seed = rng.choice([-1, 1], (m, m))
            seed = np.triu(seed)+np.triu(seed, 1).T
            for deleted in (1, 2, 3, m//2, m-1):
                key = f"m{m}_deleted{deleted}"
                observed = check(bases, seed, deleted)
                counts[key] = counts.get(key, 0)+1
                max_components[key] = max(max_components.get(key, 0), observed)
    print(json.dumps({"seed": 260906, "integer_weaves_checked": sum(counts.values()),
                      "counts": counts, "max_deleted_gram_components": max_components,
                      "status": "all reconstruction and orthogonality identities pass"}, indent=2))


if __name__ == "__main__":
    main()
