"""Exact small checks of the annealed fully-flat row count identity."""
import collections
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path
import numpy as np
from flatify_construct_2026_09_07_rank_two_weave import sylvester


records = []
for k in [2, 8]:
    h = sylvester(k)
    words = 1-2*((np.arange(1 << k, dtype=np.int64)[:, None] >> np.arange(k)) & 1)
    transformed = words @ h
    amplitude = math.isqrt(2*k)
    near = transformed[np.all((transformed == 0) | (np.abs(transformed) == amplitude), axis=1)]
    support_counts = collections.Counter(int(sum((1 << i) for i, v in enumerate(row) if v)) for row in near)
    full = (1 << k)-1
    distribution = collections.Counter()
    representatives = {}
    for perm in itertools.permutations(range(k)):
        mapped = {sum(1 << perm[i] for i in range(k) if mask & (1 << i)): count for mask, count in support_counts.items()}
        bent = sum(count*mapped.get(full ^ mask, 0) for mask, count in support_counts.items())
        distribution[bent] += 1
        representatives.setdefault(bent, perm)
    average = Fraction(sum(n*count for n, count in distribution.items()), math.factorial(k))
    predicted = Fraction(len(near)**2, math.comb(k, k//2))
    assert average == predicted
    direct = []
    parent_words = 1-2*((np.arange(1 << (2*k), dtype=np.int64)[:, None] >> np.arange(2*k)) & 1)
    for predicted_count, perm in representatives.items():
        child2 = h[list(perm)]
        parent = np.block([[h, child2], [h, -child2]])
        assert np.array_equal(parent.T @ parent, 2*k*np.eye(2*k, dtype=np.int64))
        # Inverting a permutation does not change the predicted distribution;
        # recompute its precise support transport for this row convention.
        inv = np.argsort(perm)
        mapped = {sum(1 << int(inv[i]) for i in range(k) if mask & (1 << i)): count for mask, count in support_counts.items()}
        exact_pred = sum(count*mapped.get(full ^ mask, 0) for mask, count in support_counts.items())
        value = int(np.count_nonzero(np.all(np.abs(parent_words @ parent) == math.isqrt(2*k), axis=1)))
        assert value == exact_pred
        direct.append(dict(permutation=perm, predicted=exact_pred, actual=value))
    rec = dict(k=k, near_count=len(near), support_counts=dict(support_counts),
               parent_flat_distribution=dict(distribution), average=str(average),
               predicted=str(predicted), direct_checks=direct)
    records.append(rec)
    print(json.dumps(rec), flush=True)
Path('computations/results/flatify_construct_2026_09_07_recursive_flat_count.json').write_text(
    json.dumps(dict(status='exact exhaustive permutation averages and direct Boolean checks', records=records), indent=2)+'\n')
