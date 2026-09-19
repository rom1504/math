#!/usr/bin/env python3
"""Exact local-pair rotation identities and a finite orbit-escape witness."""

import itertools
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def triangles(a):
    return np.array([a[i, j]*a[j, k]*a[k, i]
                     for i, j, k in itertools.combinations(range(len(a)), 3)], dtype=np.int64)


def rotate(a, b, r):
    cross = r[:, None] != r[None, :]
    factor = 1-2*r[:, None]*r[None, :]
    aa = np.where(cross, b, a)*factor
    bb = np.where(cross, -a, b)*factor
    return aa, bb


def permutation_works(a, b, p):
    candidate = a[np.ix_(p, p)]
    return np.array_equal(triangles(candidate), triangles(b))


def main():
    n = 4
    root_edges = list(itertools.combinations(range(1, n), 2))
    permutations = list(itertools.permutations(range(n)))
    rotations = list(itertools.product([0, 1], repeat=n))
    switches = [np.array((1,)+s, dtype=np.int64)
                for s in itertools.product([-1, 1], repeat=n-1)]
    checks = 0
    original_perm_survives = 0
    for signs in itertools.product([-1, 1], repeat=len(root_edges)):
        a = np.ones((n, n), dtype=np.int64)
        np.fill_diagonal(a, 0)
        for (i, j), sign in zip(root_edges, signs):
            a[i, j] = a[j, i] = sign
        ta = triangles(a)
        for p in permutations:
            ap = a[np.ix_(p, p)]
            for s in switches:
                b = ap*s[:, None]*s[None, :]
                tb = triangles(b)
                q = a*b
                for rr in rotations:
                    r = np.array(rr, dtype=np.int64)
                    aa, bb = rotate(a, b, r)
                    taa, tbb = triangles(aa), triangles(bb)
                    assert np.array_equal(taa*tbb, ta*tb)
                    expected_q = q*(1-2*r[:, None])*(1-2*r[None, :])
                    assert np.array_equal(aa*bb, expected_q)
                    f = ta*taa
                    original_works = permutation_works(aa, bb, p)
                    f_permuted = triangles(a[np.ix_(p, p)])*triangles(aa[np.ix_(p, p)])
                    assert original_works == np.array_equal(f, f_permuted)
                    original_perm_survives += int(original_works)
                    if np.all(ta == tb):
                        assert np.all(taa == tbb)
                    checks += 1

    a = np.array([[0,1,1,1], [1,0,-1,-1], [1,-1,0,1], [1,-1,1,0]], dtype=np.int64)
    p = (0, 2, 1, 3)
    b = a[np.ix_(p, p)]
    r = np.array([0, 0, 1, 0], dtype=np.int64)
    aa, bb = rotate(a, b, r)
    assert not any(permutation_works(aa, bb, pp) for pp in permutations)
    d = np.ones(n, dtype=np.int64)
    parent = np.block([[a, b+np.diag(d)], [b+np.diag(d), -a]])
    transform = np.eye(2*n, dtype=np.int64)
    for i in range(n):
        if r[i]:
            transform[i, i] = transform[n+i, n+i] = 0
            transform[i, n+i] = -1
            transform[n+i, i] = 1
    rotated = transform.T @ parent @ transform
    reconstructed = np.block([[aa, bb+np.diag((1-2*r)*d)],
                              [bb+np.diag((1-2*r)*d), -aa]])
    assert np.array_equal(rotated, reconstructed)
    ids = np.arange(1 << (2*n), dtype=np.int64)
    words = 1-2*((ids[:, None] >> np.arange(2*n)) & 1)
    energies = np.einsum('bi,ij,bj->b', words, parent, words)//2
    rotated_energies = np.einsum('bi,ij,bj->b', words, rotated, words)//2
    assert np.array_equal(np.sort(energies), np.sort(rotated_energies))
    result = dict(exact_identity_checks=checks,
                  original_permutation_survives=original_perm_survives,
                  example=dict(A=a.tolist(), B=b.tolist(), p=p, quarter_turns=r.tolist(),
                               A_rotated=aa.tolist(), B_rotated=bb.tolist(),
                               negative_triangle_counts=[int(np.sum(triangles(aa)<0)), int(np.sum(triangles(bb)<0))],
                               rotated_blocks_not_signed_permutation_equivalent=True,
                               all_256_parent_spin_energies_preserved=True,
                               parent_cap=int(np.max(abs(energies)))))
    target = ROOT/'computations/results/twisted_chiral_pair_rotation_audit_2026_09_19.json'
    target.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
