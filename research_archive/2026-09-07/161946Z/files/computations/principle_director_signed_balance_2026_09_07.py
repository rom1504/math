"""Exact tiny signed-incidence counts; no asymptotic inference from numerics."""
import itertools
import json
from pathlib import Path
import numpy as np


def census(n):
    edges = list(itertools.combinations(range(n), 2))
    e = len(edges)
    words = np.arange(1 << e, dtype=np.uint32)
    signs = (1 - 2 * ((words[:, None] >> np.arange(e)) & 1).astype(np.int8))
    free = [k for k, (i, j) in enumerate(edges) if i != 0]
    records = []
    for code in range(1 << len(free)):
        seed = np.ones(e, dtype=np.int8)
        for bit, k in enumerate(free):
            seed[k] = 1 - 2 * ((code >> bit) & 1)
        incidence = np.zeros((e, n), dtype=np.int8)
        for k, (i, j) in enumerate(edges):
            incidence[k, i] = 1
            incidence[k, j] = seed[k]
        sums = signs @ incidence
        count = int(np.count_nonzero(np.all(sums == 0, axis=1)))
        compatible = int(np.prod(seed)) == (-1) ** e
        assert (count > 0) == compatible
        records.append(dict(seed_code=code, parity_compatible=compatible,
                            balanced_count=count))
    positive = [r['balanced_count'] for r in records if r['balanced_count']]
    return dict(n=n, switching_classes=len(records), records=records,
                positive_min=min(positive), positive_max=max(positive))


if __name__ == '__main__':
    result = dict(status='exact finite census; parity existence verified',
                  cases=[census(3), census(5)])
    root = Path(__file__).resolve().parents[1]
    path = root / 'computations/results/principle_director_signed_balance_2026_09_07.json'
    path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k != 'cases'}))
    for case in result['cases']:
        print({k: v for k, v in case.items() if k != 'records'})
