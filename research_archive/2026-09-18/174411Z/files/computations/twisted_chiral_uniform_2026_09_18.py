#!/usr/bin/env python3
"""Exact marginal-profile selection certificates for twisted chiral doubles.

The profile certificate averages over every signed permutation, but is
computed using only two one-seed histograms. No heuristic is used in its
certificate values. Small exhaustive averages independently check the formula.
"""

import argparse
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def cube(n):
    masks = np.arange(1 << n, dtype=np.uint32)
    return (1 - 2 * ((masks[:, None] >> np.arange(n)) & 1)).astype(np.int16)


def profile_arrays(a):
    n = len(a)
    x = cube(n)
    fields = x * (x @ a)
    h = np.sum(fields, axis=1) // 2
    masks = np.arange(1 << n, dtype=np.uint32)
    cut = np.abs(h[None, :] - h[np.bitwise_xor(masks[:, None], masks[None, :])])
    internal = np.abs(x @ fields.T)
    sizes = np.sum(x == 1, axis=1)
    return h, cut, internal, sizes


def histograms(a):
    h, cut, internal, sizes = profile_arrays(a)
    n = len(a)
    length = n * (n - 1) + 1
    ch = [np.bincount(cut[sizes == k].ravel(), minlength=length) for k in range(n + 1)]
    ih = [np.bincount(internal[sizes == k].ravel(), minlength=length) for k in range(n + 1)]
    return int(np.max(np.abs(h))), ch, ih


def expected_bad(n, ch, ih, cap):
    total = Fraction(0)
    for k in range(n + 1):
        numerator = sum(
            int(ch[k][u]) * int(ih[k][v])
            for u in np.flatnonzero(ch[k])
            for v in np.flatnonzero(ih[k])
            if int(u) + int(v) > cap
        )
        total += Fraction(numerator, math.comb(n, k) * (1 << n))
    return total


def record(a, label, source):
    a = np.asarray(a, dtype=np.int16)
    n = len(a)
    q, ch, ih = histograms(a)
    curves = []
    certificate = None
    for cap in range(2 * q, 4 * n * n + 1, 2):
        bad = expected_bad(n, ch, ih, cap)
        curves.append({"cap": cap, "expected_bad": str(bad), "expected_bad_constraint_orbits": str(bad / 8)})
        if bad < 8:
            certificate = cap
            break
    target = 2 ** 1.5 * q
    rounded_target = 2 * math.floor(target / 2)
    return {
        "label": label,
        "source": source,
        "matrix": a.tolist(),
        "order": n,
        "child_cap": q,
        "core_certificate": certificate,
        "full_signing_certificate_arbitrary_diagonal": certificate + n,
        "no_error_target": target,
        "target_even_floor": rounded_target,
        "expected_bad_at_target_even_floor": str(expected_bad(n, ch, ih, rounded_target)),
        "profiles": {
            "cut": [{str(i): int(row[i]) for i in np.flatnonzero(row)} for row in ch],
            "internal": [{str(i): int(row[i]) for i in np.flatnonzero(row)} for row in ih],
        },
        "certificate_curve": curves,
    }


def exhaustive_formula_check():
    a = np.asarray([[0, 1, 1, 1], [1, 0, 1, -1], [1, 1, 0, 1], [1, -1, 1, 0]], dtype=np.int16)
    n = len(a)
    _, cut, _, _ = profile_arrays(a)
    _, ch, ih = histograms(a)
    totals = {cap: 0 for cap in range(0, 26, 2)}
    transformations = 0
    for p in itertools.permutations(range(n)):
        pa = a[np.ix_(p, p)]
        for s in cube(n):
            b = pa * s[:, None] * s[None, :]
            _, _, internal, _ = profile_arrays(b)
            value = cut + internal
            for cap in totals:
                totals[cap] += int(np.count_nonzero(value > cap))
            transformations += 1
    for cap, total in totals.items():
        assert Fraction(total, transformations) == expected_bad(n, ch, ih, cap)
    return {"order": n, "signed_permutations": transformations, "checked_thresholds": list(totals), "pass": True}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--orders", type=int, nargs="+", default=list(range(3, 11)))
    parser.add_argument("--include-chiral12", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output", default="computations/results/twisted_chiral_uniform_2026_09_18.json")
    args = parser.parse_args()
    result = {"normalization": "Q(A)=max |sum_(i<j) A_ij x_i x_j|", "records": []}
    if args.check:
        result["independent_exhaustive_check"] = exhaustive_formula_check()
    for n in args.orders:
        source = ROOT / f"computations/results/m{n}_minimizer_orbits.json"
        if source.exists():
            data = json.loads(source.read_text())
            seeds = [(f"m{n}_class{item['class']}", item["representative_matrix"]) for item in data["classes"]]
        else:
            source = ROOT / f"computations/results/exact_m{n}.json"
            data = json.loads(source.read_text())
            seeds = [(f"m{n}_stored", data["matrix"])]
        for label, a in seeds:
            item = record(a, label, str(source.relative_to(ROOT)))
            result["records"].append(item)
            print(label, item["child_cap"], item["core_certificate"], item["expected_bad_at_target_even_floor"], flush=True)
    if args.include_chiral12:
        source = ROOT / "artifacts/dependent_profile_recovery_witness.json"
        data = json.loads(source.read_text())
        item = record(data["matrix"], "compressed_chiral12", str(source.relative_to(ROOT)))
        result["records"].append(item)
        print(item["label"], item["child_cap"], item["core_certificate"], item["expected_bad_at_target_even_floor"], flush=True)
    output = ROOT / args.output
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
