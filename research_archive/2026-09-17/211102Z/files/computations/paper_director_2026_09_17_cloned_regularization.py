"""Exact finite check of amplified fresh-block and nearcode-width identities.

Enumerates every random block signing and every physical spin. All
inequalities tested below use integers, not floating-point tolerance.
This checks the construction, not the asymptotic probability proof.
"""
import itertools
import json
import numpy as np


def cube(d):
    return np.array(list(itertools.product((-1, 1), repeat=d)), dtype=np.int64)


def keys(m, p):
    return [("old", i, a) for a in range(p) for i in range(m)] + [
        ("block", a, b) for a in range(p) for b in range(a + 1, p)]


def matrix(m, p, r, labels, signs):
    n = m + p*r
    a = np.ones((n, n), dtype=np.int64) - np.eye(n, dtype=np.int64)
    # A nonconstant old core when available; same fixed core in every world.
    for i in range(m):
        for j in range(i+1, m):
            a[i, j] = a[j, i] = -1 if (i+j) % 3 == 0 else 1
    for label, sign in zip(labels, signs):
        kind, i, j = label
        if kind == "old":
            a[i, m+j*r:m+(j+1)*r] = sign
            a[m+j*r:m+(j+1)*r, i] = sign
        else:
            a[m+i*r:m+(i+1)*r, m+j*r:m+(j+1)*r] = sign
            a[m+j*r:m+(j+1)*r, m+i*r:m+(i+1)*r] = sign
    return a


def run(m, p, r):
    base_keys, ext_keys = keys(m, p), keys(m, p+1)
    extract = [ext_keys.index(key) for key in base_keys]
    base_spins, ext_spins, rows = cube(m+p*r), cube(m+(p+1)*r), cube(m)
    conditional_sums = {}
    conditional_counts = {}
    for signs in cube(len(ext_keys)):
        a = matrix(m, p+1, r, ext_keys, signs)
        energies = np.einsum("bi,ij,bj->b", ext_spins, a, ext_spins)//2
        base = tuple(int(signs[i]) for i in extract)
        conditional_sums[base] = conditional_sums.get(base, 0) + int(np.abs(energies).max())
        conditional_counts[base] = conditional_counts.get(base, 0) + 1
    checked = 0
    minimum_margin = None
    for signs in cube(len(base_keys)):
        a = matrix(m, p, r, base_keys, signs)
        energies = np.einsum("bi,ij,bj->b", base_spins, a, base_spins)//2
        cap = int(np.abs(energies).max())
        base = tuple(int(s) for s in signs)
        count = conditional_counts[base]
        delta_numerator = conditional_sums[base] - count*cap
        assert delta_numerator >= 0
        for tolerance in sorted({0, 1, 2, 4, cap}):
            near = base_spins[cap-np.abs(energies) <= tolerance]
            projected = np.unique(near[:, :m], axis=0)
            width_numerator = int(np.max(rows @ projected.T, axis=1).sum())
            left = r*width_numerator*count
            right = ((tolerance+r*(r-1)//2)*count+delta_numerator)*len(rows)
            assert left <= right, (m,p,r,base,tolerance,left,right)
            margin = right-left
            minimum_margin = margin if minimum_margin is None else min(minimum_margin,margin)
            checked += 1
    return {"m":m,"p":p,"r":r,"base_random_worlds":2**len(base_keys),
            "extended_random_worlds":2**len(ext_keys),
            "physical_spin_count":len(base_spins),"all_window_checks":checked,
            "minimum_integer_margin":minimum_margin}


if __name__ == "__main__":
    cases = [(2,1,1),(2,1,2),(2,2,2),(3,1,2),(2,1,3),(2,2,3)]
    reports = []
    for case in cases:
        result = run(*case)
        reports.append(result)
        print(json.dumps(result), flush=True)
    print(json.dumps({"status":"PASS exact amplified fresh-block identity", "cases":reports},indent=2))
