#!/usr/bin/env python3
"""Independent finite audit of the paired-polar inequality.

Enumerated energies, bilinear caps, and cut weights are exact integers.
Eigenvalues/arcsines are floating-point diagnostics, not proof certificates.
The asymptotic/theorem claims rely on the accompanying analytic proof.
No files are written. Python 3.9 compatible.
"""
import argparse
import itertools
import json
import math

import numpy as np


def spin_table(n):
    masks = np.arange(1 << (n - 1), dtype=np.uint64)
    shifts = np.arange(n, dtype=np.uint64)
    return (1 - 2 * ((masks[:, None] >> shifts) & 1).astype(np.int64))


def audit_order(n, batch_size=256):
    edges = list(itertools.combinations(range(1, n), 2))
    spins = spin_table(n)
    minimum = n * (n - 1) // 2
    minimum_classes = 0
    worst_slack = float('inf')
    worst_nuclear_half_deficit = -float('inf')
    total = 1 << len(edges)
    for start in range(0, total, batch_size):
        masks = np.arange(start, min(total, start + batch_size), dtype=np.uint64)
        matrices = np.ones((len(masks), n, n), dtype=np.int64)
        for i in range(n):
            matrices[:, i, i] = 0
        for bit, (i, j) in enumerate(edges):
            value = 1 - 2 * ((masks >> bit) & 1).astype(np.int64)
            matrices[:, i, j] = value
            matrices[:, j, i] = value
        fields = np.einsum('bij,sj->bsi', matrices, spins)
        energy_twice = np.einsum('si,bsi->bs', spins, fields)
        assert np.all(energy_twice % 2 == 0)
        caps = np.max(np.abs(energy_twice), axis=1) // 2
        bilinear = np.max(np.sum(np.abs(fields), axis=2), axis=1)
        assert np.all(2 * caps <= bilinear)
        assert np.all(bilinear <= 4 * caps)
        eigs = np.linalg.eigvalsh(matrices.astype(float))
        nuclear = np.sum(np.abs(eigs), axis=1)
        ratios = n / nuclear
        assert np.max(ratios) <= 1 + 1e-10
        polar_lower = n * (n - 1) / math.pi * np.arcsin(np.minimum(ratios, 1))
        slack = caps - polar_lower
        assert np.min(slack) >= -1e-9
        worst_slack = min(worst_slack, float(np.min(slack)))
        worst_nuclear_half_deficit = max(
            worst_nuclear_half_deficit,
            float(np.max((nuclear / 2 - caps) / n)))
        found = int(np.min(caps))
        if found < minimum:
            minimum = found
            minimum_classes = int(np.sum(caps == found))
        elif found == minimum:
            minimum_classes += int(np.sum(caps == found))
    return dict(n=n, root_gauged_signings=total, exact_M=minimum,
                root_gauged_minimizers=minimum_classes,
                least_polar_slack_numeric=worst_slack,
                largest_nuclear_half_deficit_per_vertex_numeric=
                worst_nuclear_half_deficit)


def cut_laplacian_counterexample():
    a = np.array([[0, 1, 1, 1, 1], [1, 0, 1, -1, 1],
                  [1, 1, 0, 1, -1], [1, -1, 1, 0, -1],
                  [1, 1, -1, -1, 0]], dtype=np.int64)
    x = spin_table(5)
    energies = np.einsum('si,ij,sj->s', x, a, x) // 2
    cap = int(np.max(np.abs(energies)))
    row = a.sum(axis=1)
    lap = np.diag(row) - a
    cuts = (int(a.sum() // 2) - energies) // 2
    assert cap == 4 and int(np.min(cuts)) == 0 and int(np.max(cuts)) == 4
    det = int(lap[0, 0] * lap[3, 3] - lap[0, 3] ** 2)
    assert det == -1
    return dict(matrix=a.tolist(), exact_cap=cap, cut_range=[0, 4],
                negative_principal_minor_indices=[0, 3],
                exact_minor_determinant=det)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-order', type=int, default=7)
    args = parser.parse_args()
    assert 2 <= args.max_order <= 8
    result = dict(
        classification='exact exhaustive caps; floating-point analytic-bound diagnostics',
        polar_statement='Q(A) >= n(n-1)/pi * asin(n / nuclear_norm(A))',
        enumerations=[audit_order(n) for n in range(2, args.max_order + 1)],
        exact_cut_laplacian_counterexample=cut_laplacian_counterexample())
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
