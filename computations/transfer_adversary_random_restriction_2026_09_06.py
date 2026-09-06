"""Finite exact checks for the random-restriction theorem.

The asymptotic assertion is proved in the companion artifact; no numerical
extrapolation is used as its evidence.
"""

from collections import Counter
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def exact_greedy_increment_law(k):
    edges = [(i, j) for i in range(1, k) for j in range(i)]
    counts = Counter()
    for bits in range(1 << len(edges)):
        a = np.zeros((k, k), dtype=np.int64)
        for bit, (i, j) in enumerate(edges):
            a[i, j] = a[j, i] = 1 - 2 * ((bits >> bit) & 1)
        x = np.ones(k, dtype=np.int64)
        zs = []
        for i in range(1, k):
            field = int(a[i, :i] @ x[:i])
            x[i] = 1 if field >= 0 else -1
            zs.append(abs(field))
        assert int(x @ a @ x) // 2 == sum(zs)
        counts[tuple(zs)] += 1

    independent_counts = [Counter() for _ in range(k - 1)]
    for length in range(1, k):
        for plus in range(length + 1):
            independent_counts[length - 1][abs(2 * plus - length)] += math.comb(length, plus)
    for zs, count in counts.items():
        target = math.prod(independent_counts[i][z] for i, z in enumerate(zs))
        assert count == target
    assert len(counts) == math.prod(map(len, independent_counts))
    mean = Fraction(sum(sum(zs) * count for zs, count in counts.items()), 1 << len(edges))
    return {"k": k, "matrices_enumerated": 1 << len(edges),
            "increment_joint_atoms": len(counts), "independence_exact": True,
            "greedy_mean": str(mean)}


def finite_restriction_check():
    a = np.array(json.loads((ROOT / "computations/results/exact_m8.json").read_text())["matrix"])
    n = len(a)
    k = 3
    edges = list(itertools.combinations(range(k), 2))
    counts = Counter()
    for labels in itertools.permutations(range(n), k):
        code = sum(int(a[labels[i], labels[j]] == 1) << bit
                   for bit, (i, j) in enumerate(edges))
        counts[code] += 1
    total = math.perm(n, k)
    tv = sum(abs(Fraction(counts[code], total) - Fraction(1, 1 << len(edges)))
             for code in range(1 << len(edges))) / 2
    indicators = np.array(list(itertools.product((0, 1), repeat=n)), dtype=np.int64)
    cut_integer = int(np.max(np.abs(indicators @ a @ indicators.T)))
    delta = Fraction(cut_integer, 2 * n * n)
    bound = (1 << (len(edges) - 1)) * len(edges) * delta + Fraction(len(edges), n)
    assert tv <= bound
    return {"seed_order": n, "restriction_order": k,
            "exact_tv_from_iid": str(tv), "exact_graphon_cut_norm": str(delta),
            "finite_tv_upper_bound": str(bound)}


def main():
    print(json.dumps({"iid_check": exact_greedy_increment_law(6),
                      "restriction_check": finite_restriction_check(),
                      "greedy_asymptotic_constant": (2 / 3) * math.sqrt(2 / math.pi)},
                     indent=2))


if __name__ == "__main__":
    main()
