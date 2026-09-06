"""Finite existential-partition screens, frozen before inspecting results.

BP-half: at N>=10, some half partition has both normalized child caps<1/2.
BP-exact: some half partition has both children exactly minimizing.
RF-partition-1: for each prescribed comparable m,n, some partition obeys
  Q(child1)^(2/3)+Q(child2)^(2/3)<=Q(parent)^(2/3)+sqrt(N).

Only stored representatives are tested. Their global minimality status is
imported from the provenance file, not reproved here. Every cap and each
inequality comparison is replayed with integers or rational intervals.
No output files are written.
"""

import itertools
import json
import argparse
from collections import Counter
from functools import lru_cache
from math import isqrt
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
MINIMUM = (0, 0, 1, 3, 4, 4, 5, 9, 10, 12, 13, 17, 18, 20, 21)
SCALE = 10**9


@lru_cache(None)
def spins(n):
    return np.array([(1,) + x for x in itertools.product((-1, 1), repeat=n - 1)], dtype=np.int64)


def cap(a):
    x = spins(len(a))
    twice = np.einsum("bi,ij,bj->b", x, a, x)
    assert not np.any(twice % 2)
    return int(np.max(np.abs(twice))) // 2


@lru_cache(None)
def powered_interval(q):
    target = q * q * SCALE**3
    lo, hi = 0, (q + 1) * SCALE
    while hi - lo > 1:
        middle = (lo + hi) // 2
        if middle**3 <= target:
            lo = middle
        else:
            hi = middle
    return lo, lo if lo**3 == target else lo + 1


metadata = json.loads((ROOT / "computations/results/transfer_adversary_minimizer_isotropy_2026_09_06.json").read_text())
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--external-witnesses", action="store_true",
                    help="Also test Q15=27 and Q16=30 witnesses, WITHOUT claiming their exact optimality.")
options = parser.parse_args()
if options.external_witnesses:
    path = "computations/results/transfer_adversary_external_order15_16_witness_verify_2026_09_06.json"
    for index, (n, q) in enumerate(((15, 27), (16, 30))):
        metadata["cases"].append({
            "case": f"n{n}_external_witness_NOT_optimality_replayed",
            "source": {"path": path, "json_pointer": f"/cases/{index}/matrix",
                       "minimum_cap_imported": q},
        })
total_partitions = 0
for item in metadata["cases"]:
    provenance = item["source"]
    matrix = json.loads((ROOT / provenance["path"]).read_text())
    for piece in provenance["json_pointer"].strip("/").split("/"):
        matrix = matrix[int(piece)] if isinstance(matrix, list) else matrix[piece]
    a = np.array(matrix, dtype=np.int64)
    n = len(a)
    q = cap(a)
    assert q == provenance["minimum_cap_imported"]
    if n < len(MINIMUM):
        assert q == MINIMUM[n]
    balanced = None
    rf_results = []
    for m in range((n + 2) // 3, n // 2 + 1):
        pairs = Counter()
        witnesses = {}
        for left in itertools.combinations(range(n), m):
            if m * 2 == n and 0 not in left:
                continue
            right = tuple(i for i in range(n) if i not in left)
            pair = cap(a[np.ix_(left, left)]), cap(a[np.ix_(right, right)])
            pairs[pair] += 1
            witnesses.setdefault(pair, left)
            total_partitions += 1
        half_count = sum(count for (u, v), count in pairs.items()
                         if 4 * u * u < m**3 and 4 * v * v < (n - m)**3)
        exact_count = pairs[MINIMUM[m], MINIMUM[n - m]]
        lower_parent, upper_parent = powered_interval(q)
        lower_sqrt = isqrt(n * SCALE**2)
        upper_sqrt = lower_sqrt if lower_sqrt**2 == n * SCALE**2 else lower_sqrt + 1
        passes = [pair for pair in pairs if powered_interval(pair[0])[1]
                  + powered_interval(pair[1])[1] - lower_parent <= lower_sqrt]
        definitely_fails = all(powered_interval(u)[0] + powered_interval(v)[0]
                               - upper_parent > upper_sqrt for u, v in pairs)
        rf_results.append((m, "PASS" if passes else "FAIL" if definitely_fails else "UNRESOLVED_INTERVAL",
                           "both_exact_count", exact_count))
        if m == n // 2:
            balanced = {"partitions": sum(pairs.values()), "strict_subhalf_count": half_count,
                        "both_exact_count": exact_count, "pairs": dict(sorted(pairs.items()))}
    print(item["case"], "cap", q, "balanced", balanced, "RF1", rf_results)

print("Total exact partition evaluations:", total_partitions)
print("Finite representative screen only; no asymptotic or all-minimizer conclusion.")
