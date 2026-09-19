#!/usr/bin/env python3
"""Complete small-order average of the stable-permanent selector.

Exact rational Brégman upper bounds use certified integer ceilings of
factorial roots. This is exhaustive for the chosen n3/n4 children, not an
asymptotic or global signing-optimality certificate.
"""
import itertools
import json
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path

import numpy as np

from twisted_chiral_stable_permanent_2026_09_19 import convolve, polynomial_permanent

ROOT = Path(__file__).resolve().parents[1]
SCALE = 10 ** 6


def factorial_root_ceiling(r):
    if r == 0:
        return 0
    target = math.factorial(r) * SCALE ** r
    lo, hi = 0, r * SCALE
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if mid ** r >= target:
            hi = mid
        else:
            lo = mid
    assert hi ** r >= target and (hi - 1) ** r < target
    return hi


def fraction_record(num, den):
    q = Fraction(num, den)
    return {"numerator": q.numerator, "denominator": q.denominator, "value": float(q)}


def run(n):
    A = np.array(json.loads((ROOT / ("computations/results/exact_m%d.json" % n)).read_text())["matrix"], dtype=np.int64)
    spins = np.array(list(itertools.product([-1, 1], repeat=n)), dtype=np.int64)
    energies = np.einsum("bi,ij,bj->b", spins, A, spins) // 2
    fields = spins @ A
    roots = [factorial_root_ceiling(i) for i in range(n + 1)]
    raw, stable, upper = Counter(), Counter(), Counter()
    cases = 0
    for ix, x in enumerate(spins):
        if x[0] != 1:
            continue
        alpha = x * fields[ix]
        for iy, y in enumerate(spins):
            gamma = y * fields[iy]
            plus = np.flatnonzero(x * y == 1).tolist()
            minus = np.flatnonzero(x * y == -1).tolist()
            count_p = math.factorial(len(plus)) * math.factorial(len(minus))
            for z in spins:
                for chosen in itertools.combinations(range(n), len(plus)):
                    image_plus = list(chosen)
                    image_minus = [j for j in range(n) if j not in chosen]
                    relative = -np.ones(n, dtype=np.int64)
                    relative[image_plus] = 1
                    w = z * relative
                    b, c = z * (A @ w), w * (A @ z)
                    margins = np.minimum(alpha[:, None] + b[None, :],
                                         -gamma[:, None] + c[None, :])
                    numerator = convolve(polynomial_permanent(plus, image_plus, margins),
                                         polynomial_permanent(minus, image_minus, margins))
                    base = int(energies[ix] - energies[iy] + z @ A @ w)
                    for v, value in numerator.items():
                        stable[base + v] += value
                    bounded = Counter({0: 1})
                    for rows, cols in [(plus, image_plus), (minus, image_minus)]:
                        for i in rows:
                            d0 = sum(int(margins[i, j] >= 0) for j in cols)
                            d2 = sum(int(margins[i, j] >= 2) for j in cols)
                            bounded = convolve(bounded, {1: roots[d0], -1: roots[d2]})
                    for v, value in bounded.items():
                        upper[base + v] += value
                    for positive in range(n + 1):
                        raw[base + 2 * positive - n] += count_p * math.comb(n, positive)
                    cases += 1
    denominator = (4 ** n) * math.factorial(n)
    assert sum(raw.values()) == denominator * 2 ** (2 * n - 1)
    assert all(upper[e] >= value * SCALE ** n for e, value in stable.items())
    # Independent full twist/matching enumeration: no conditional fields,
    # sector permanents, or row-degree polynomial is used in this audit.
    X = np.array([(1,) + tail for tail in itertools.product([-1, 1], repeat=2*n-1)], dtype=np.int64)
    direct_raw, direct_stable = Counter(), Counter()
    direct_parents = 0
    for p_tuple in itertools.permutations(range(n)):
        p = np.array(p_tuple, dtype=np.int64)
        for s in spins:
            B = A[np.ix_(p, p)] * np.outer(s, s)
            for d in spins:
                C = B + np.diag(d)
                D = np.block([[A, C], [C, -A]])
                fields_direct = (X @ D) * X
                e = fields_direct.sum(axis=1) // 2
                local = np.all(fields_direct > 0, axis=1)
                for value, count in zip(*np.unique(e, return_counts=True)):
                    direct_raw[int(value)] += int(count)
                for value, count in zip(*np.unique(e[local], return_counts=True)):
                    direct_stable[int(value)] += int(count)
                direct_parents += 1
    assert direct_parents == denominator
    assert direct_raw == raw and direct_stable == stable
    thresholds = []
    for L in range(0, n * (2 * n - 1) + 1):
        sr = sum(v for e, v in stable.items() if e > L)
        rr = sum(v for e, v in raw.items() if e > L)
        ur = sum(v for e, v in upper.items() if e > L)
        thresholds.append({"L": L,
                           "exact_positive_raw_above": fraction_record(rr, denominator),
                           "exact_positive_stable_above": fraction_record(sr, denominator),
                           "certified_bregman_upper_above": fraction_record(ur, denominator * SCALE ** n),
                           "exact_stable_selector_certifies": sr < denominator,
                           "bregman_selector_certifies": ur < denominator * SCALE ** n})
    return {"n": n, "matrix": A.tolist(), "child_cap": int(np.max(np.abs(energies))),
            "conditional_cases": cases, "group_matching_denominator": denominator,
            "independent_direct_parent_count": direct_parents,
            "independent_full_parent_histogram_audit_pass": True,
            "factorial_root_upper_scale": SCALE, "factorial_root_upper_numerators": roots,
            "positive_raw_energy_histogram_numerators": dict(sorted(raw.items())),
            "positive_stable_energy_histogram_numerators": dict(sorted(stable.items())),
            "bregman_energy_histogram_upper_numerators": dict(sorted(upper.items())),
            "thresholds": thresholds}


def main():
    records = [run(n) for n in (3, 4)]
    output = {"status": "complete exact finite average, with rational upper certificates",
              "positive_state_convention": "projective: x0=+1, all y; criterion cutoff1",
              "records": records}
    path = ROOT / "computations/results/twisted_chiral_stable_average_2026_09_19.json"
    path.write_text(json.dumps(output, indent=2) + "\n")
    for rec in records:
        fields = ["exact_positive_raw_above", "exact_positive_stable_above", "certified_bregman_upper_above"]
        print(json.dumps({"n": rec["n"], "cases": rec["conditional_cases"],
                          "first_L_with_count_below1": {field: next(row["L"] for row in rec["thresholds"] if row[field]["numerator"] < row[field]["denominator"]) for field in fields}}))


if __name__ == "__main__":
    main()
