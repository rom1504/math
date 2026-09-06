"""Exact small balanced-selector couplings and weave cap continuity.

No input/output files. Uses only finite integer/rational checks and numpy.
"""
from collections import defaultdict
from fractions import Fraction as F
import itertools
import json
import math

import numpy as np


def repair_choices(selected, fibre, k):
    if len(selected) >= k:
        return [frozenset(x) for x in itertools.combinations(selected, k)]
    missing = sorted(set(fibre) - set(selected))
    return [frozenset(selected) | frozenset(x)
            for x in itertools.combinations(missing, k - len(selected))]


def exact_coupling(m, k):
    size, retained = m * m, m * k
    fibres = [tuple(range(i * m, (i + 1) * m)) for i in range(m)]
    law = defaultdict(F)
    mean_h = F(0)
    selections = math.comb(size, retained)
    for raw in itertools.combinations(range(size), retained):
        selected = frozenset(raw)
        choices = [repair_choices(selected.intersection(fibre), fibre, k)
                   for fibre in fibres]
        h = sum(abs(len(selected.intersection(fibre)) - k) for fibre in fibres)
        mean_h += F(h, selections)
        weight = F(1, selections * math.prod(map(len, choices)))
        for pieces in itertools.product(*choices):
            repaired = frozenset().union(*pieces)
            assert len(selected.symmetric_difference(repaired)) == h
            law[repaired] += weight
    balanced = math.comb(m, k) ** m
    assert len(law) == balanced
    assert set(law.values()) == {F(1, balanced)}
    variance = F(k) * F(m - 1, m) * F(size - retained, size - 1)
    assert mean_h ** 2 <= m * m * variance
    return {'m': m, 'k': k, 'uniform_selectors': selections,
            'balanced_selectors': balanced, 'mean_symmetric_difference': str(mean_h)}


def weave(m, rng):
    base = np.array([[1 - 2 * (bin(i & j).count('1') % 2)
                      for j in range(m)] for i in range(m)], dtype=np.int64)
    bases = []
    for _ in range(m):
        h = base[rng.permutation(m)][:, rng.permutation(m)]
        h = h * rng.choice((-1, 1), size=(m, 1)) * rng.choice((-1, 1), size=(1, m))
        bases.append(h)
    upper = np.triu(rng.choice((-1, 1), size=(m, m)))
    signs = upper + upper.T - np.diag(np.diag(upper))
    out = np.empty((m * m, m * m), dtype=np.int64)
    for i, j, a, b in itertools.product(range(m), repeat=4):
        out[i * m + a, j * m + b] = signs[i, j] * bases[i][a, j] * bases[j][b, i]
    assert np.array_equal(out, out.T)
    assert np.array_equal(out @ out, m * m * np.eye(m * m, dtype=np.int64))
    np.fill_diagonal(out, 0)
    return out


def raw_cap(a):
    n = len(a)
    spins = np.array([(1,) + x for x in itertools.product((-1, 1), repeat=n - 1)], dtype=np.int64)
    doubled = int(np.max(np.abs(np.einsum('bi,ij,bj->b', spins, a, spins))))
    assert doubled % 2 == 0
    return doubled // 2


def cap_checks():
    rng = np.random.default_rng(20260906)
    count = 0
    for m, k in ((2, 1), (4, 1), (4, 2), (4, 3)):
        for _ in range(25):
            a = weave(m, rng)
            selected = frozenset(map(int, rng.choice(m * m, size=m * k, replace=False)))
            pieces = []
            for i in range(m):
                fibre = tuple(range(i * m, (i + 1) * m))
                choices = repair_choices(selected.intersection(fibre), fibre, k)
                pieces.append(choices[int(rng.integers(len(choices)))])
            repaired = frozenset().union(*pieces)
            h = len(selected.symmetric_difference(repaired))
            left, right = sorted(selected), sorted(repaired)
            gap = abs(raw_cap(a[np.ix_(left, left)]) - raw_cap(a[np.ix_(right, right)]))
            assert gap * gap <= (m + 1) ** 2 * (m * k) * h
            count += 1
    return count


def main():
    print(json.dumps({'exact_couplings': [exact_coupling(m, k) for m, k in
                                         ((2, 1), (3, 1), (3, 2), (4, 1))],
                      'integer_nearby_selector_cap_checks': cap_checks(),
                      'seed': 20260906}, indent=2))


if __name__ == '__main__':
    main()
