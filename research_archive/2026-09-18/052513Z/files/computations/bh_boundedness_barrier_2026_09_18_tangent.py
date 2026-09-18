#!/usr/bin/env python3
"""Independent exact audit of the localization agent's degree-transfer tangent.

All coefficients and pointwise identities are checked in Z[zeta_6].
Floating-point BH ratios at the end are diagnostics only.
"""
from collections import Counter
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path


def add(x, y):
    return x[0]+y[0], x[1]+y[1]


def neg(x):
    return -x[0], -x[1]


def mul(x, y):
    a, b = x
    c, d = y
    return a*c-b*d, a*d+b*c+b*d


def conj(x):
    return x[0]+x[1], -x[1]


def norm2(x):
    return x[0]**2+x[0]*x[1]+x[1]**2


def degree(mask):
    return bin(mask).count("1")


def convolution(x, y):
    result = {}
    for i, a in x.items():
        for j, b in y.items():
            result[i ^ j] = add(result.get(i ^ j, (0, 0)), mul(a, b))
    return {i: a for i, a in result.items() if a != (0, 0)}


def evaluate(p, x):
    out = (0, 0)
    for mask, value in p.items():
        out = add(out, neg(value) if degree(mask & x) % 2 else value)
    return out


def tensor(x, y):
    return {i+32*j: mul(a, b) for i, a in x.items() for j, b in y.items()}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    roots = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
    masks = [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 16, 17, 18, 20, 24]
    exponents = [0, 2, 0, 5, 1, 4, 2, 4, 3, 1, 4, 1, 2, 2, 3, 2]
    fnum = {i: roots[e] for i, e in zip(masks, exponents)}  # denominator 4
    anum = {0: conj(roots[1]), 4: roots[1], 8: (1, 0)}  # denominator 1
    bnum = convolution(convolution(fnum, fnum),
                       {i: conj(a) for i, a in anum.items()})  # denominator 16
    assert max(map(degree, bnum)) == 3
    z_histogram = Counter()
    for x in range(32):
        f, a, b = evaluate(fnum, x), evaluate(anum, x), evaluate(bnum, x)
        assert f in [(4*c, 4*d) for c, d in roots]
        assert b == mul(mul(f, f), conj(a))
        znum = mul(a, conj(f))  # z denominator 4
        assert znum in [(0, 0)] + [(8*c, 8*d) for c, d in roots]
        z_histogram[znum] += 1
    Fnum = tensor(fnum, fnum)  # denominator 16
    hnum = tensor(anum, bnum)
    for mask, value in tensor(bnum, anum).items():
        hnum[mask] = add(hnum.get(mask, (0, 0)), neg(value))
    hnum = {i: value for i, value in hnum.items() if value != (0, 0)}
    assert max(map(degree, hnum)) <= 4
    assert hnum[32*21] == (12, -12)
    assert 32*21 not in Fnum
    inner = (0, 0)
    for mask, value in Fnum.items():
        inner = add(inner, mul(conj(value), hnum.get(mask, (0, 0))))
    assert 2*inner[0]+inner[1] == 0
    ratio_histogram = Counter()
    for x in range(1024):
        F, h = evaluate(Fnum, x), evaluate(hnum, x)
        assert norm2(F) == 256
        ratio_num = mul(h, conj(F))
        assert ratio_num in [(0, 0), (-1024, 2048), (1024, -2048)]
        ratio_histogram[ratio_num] += 1
    outside = {mask: value for mask, value in hnum.items() if mask not in Fnum}
    outside_norm_counts = Counter(norm2(value) for value in outside.values())
    assert outside_norm_counts == {48: 24, 144: 12}
    assert sum(norm2(value) for value in outside.values()) == 2880
    h2num = convolution(hnum, {i: conj(a) for i, a in hnum.items()})
    assert max(map(degree, h2num)) == 8  # denominator 256
    assert convolution(h2num, h2num) == {
        i: (12288*a[0], 12288*a[1]) for i, a in h2num.items()}
    assert convolution(hnum, h2num) == {
        i: (12288*a[0], 12288*a[1]) for i, a in hnum.items()}
    Cnum = convolution(Fnum, h2num)  # C=F|h|^2, denominator 4096
    assert max(map(degree, Cnum)) == 8
    assert Cnum[503] == (3072, 0)
    assert sum(map(norm2, hnum.values())) == 18*256
    assert sum(map(norm2, Cnum.values())) == 864*4096**2
    assert 7*Fraction(81, 4)+256*Fraction(153, 256)**2 < 256
    delta = Fraction(1, 2**16)
    quotient_excess = (Fraction(45, 2**18)-4*Fraction(3, 81920))/(1+Fraction(3, 81920))
    assert quotient_excess == Fraction(33, 16*81923)
    assert quotient_excess > delta
    assert 31**5 > 256*10**5
    assert Fraction(31, 10)+delta < 4
    q = 8/5
    outside_q_mass = sum((math.sqrt(norm2(value))/16)**q
                         for value in outside.values())
    diagnostics = []
    for denominator in (10000, 1000, 100, 50, 25, 20, 10):
        all_masks = set(Fnum) | set(hnum)
        q_mass = 0.0
        for mask in all_masks:
            a = Fnum.get(mask, (0, 0))
            b = hnum.get(mask, (0, 0))
            numerator = add((denominator*a[0], denominator*a[1]), b)
            q_mass += (math.sqrt(norm2(numerator))/(16*denominator))**q
        cap = math.sqrt(1+48/denominator**2)
        ratio = q_mass**(1/q)/cap
        diagnostics.append({"epsilon": 1/denominator, "ratio": ratio,
                            "convexity_lower_ratio":
                                (2**q+outside_q_mass/denominator**q)**(1/q)/cap})
    result = {"exact_status": "PASS", "base_degree": 2,
                      "a_degree": 1, "b_degree": 3,
                      "F_degree": 4, "h_degree": 4,
                      "F_support": len(Fnum), "h_support": len(hnum),
                      "new_support": len(outside),
                      "h_squared_degree": 8,
                      "C_F_times_h_squared_degree": 8,
                      "C_coefficient_mask503": "3/4",
                      "normalization_idempotence_and_absorption": "PASS",
                      "h_squared_L2": 18,
                      "C_squared_L2": 864,
                      "outside_coefficient_squared_numerator_counts_denominator256":
                          dict(outside_norm_counts),
                      "outside_squared_mass": "45/4",
                      "certified_degree4_ratio_lower_bound": "2 + 2^(-18), strictly",
                      "b_denominator": 16,
                      "b_numerators": {str(k): list(v) for k, v in sorted(bnum.items())},
                      "z_histogram_denominator4": {str(k): v for k, v in z_histogram.items()},
                      "tangent_ratio_histogram_denominator256":
                          {str(k): v for k, v in ratio_histogram.items()},
                      "outside_q_mass_diagnostic": outside_q_mass,
                      "ratio_diagnostics_not_certificates": diagnostics}
    payload = json.dumps(result, indent=2) + "\n"
    if args.output:
        args.output.write_text(payload)
    print(payload)


if __name__ == "__main__":
    main()
