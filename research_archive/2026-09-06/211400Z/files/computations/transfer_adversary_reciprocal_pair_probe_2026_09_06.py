"""Exact small Walsh reciprocal-spectrum probe, not an asymptotic proof."""

import itertools
import json
from fractions import Fraction
from pathlib import Path

import numpy as np


def hadamard(n):
    h = np.ones((1, 1), dtype=np.int64)
    while len(h) < n:
        h = np.block([[h, h], [h, -h]])
    return h


def probe(n):
    h = hadamard(n)
    codes = np.arange(1 << (n - 1), dtype=np.uint64)
    tail = (1 - 2 * ((codes[:, None] >> np.arange(n - 1, dtype=np.uint64)) & 1)).astype(np.int64)
    fs = np.column_stack((np.ones(len(codes), dtype=np.int64), tail))
    af = fs @ h
    square = af * af
    number_flips = n // 8
    best = 10 ** 30
    best_data = None
    best_l1 = 10 ** 30
    best_l1_data = None
    for positions in itertools.combinations(range(n), number_flips):
        bf = af.copy()
        for position in positions:
            bf -= 2 * fs[:, position, None] * h[position]
        numerators = np.sum(square * bf * bf, axis=1)
        product_target = n - 2 * number_flips
        l1_numerators = np.sum(np.abs(af * bf - product_target), axis=1)
        l1_idx = int(np.argmin(l1_numerators))
        if int(l1_numerators[l1_idx]) < best_l1:
            best_l1 = int(l1_numerators[l1_idx])
            f = fs[l1_idx].copy()
            g = f.copy()
            g[list(positions)] *= -1
            best_l1_data = {"f": f.tolist(), "g": g.tolist(),
                            "walsh_f": af[l1_idx].tolist(), "walsh_g": bf[l1_idx].tolist()}
        idx = int(np.argmin(numerators))
        value = int(numerators[idx])
        if value < best:
            best = value
            f = fs[idx].copy()
            g = f.copy()
            g[list(positions)] *= -1
            best_data = {"f": f.tolist(), "g": g.tolist(),
                         "walsh_f": af[idx].tolist(), "walsh_g": bf[idx].tolist()}
    rho = Fraction(n - 2 * number_flips, n)
    return {"group_order": n, "correlation": str(rho),
            "pairs_examined_up_to_common_reversal": len(codes) * __import__("math").comb(n, number_flips),
            "min_fourier_product_second_moment": str(Fraction(best, n ** 3)),
            "cauchy_lower_target": str(rho * rho),
            "min_reciprocal_product_l1_defect": str(Fraction(best_l1, n ** 2)),
            "below_absolute_correlation": Fraction(best, n ** 3) < abs(rho),
            "l1_witness": best_l1_data,
            "witness": best_data}


if __name__ == "__main__":
    witness_path = Path(__file__).with_name("results") / "transfer_adversary_reciprocal_order64_witness_2026_09_06.json"
    with witness_path.open() as witness_file:
        witness = json.load(witness_file)
    f, g = np.array(witness["f"], dtype=np.int64), np.array(witness["g"], dtype=np.int64)
    h64 = hadamard(64)
    assert int(f @ g) == witness["input_inner_product"] == 48
    assert np.array_equal(h64 @ f, witness["walsh_f"])
    assert np.array_equal(h64 @ g, witness["walsh_g"])
    numerator = int(np.sum(np.abs((h64 @ f) * (h64 @ g) - 48)))
    assert numerator == witness["integer_walsh_product_l1_numerator"] == 1312
    print(json.dumps({"exhaustive": [probe(8), probe(16)],
                      "fixed_order64_witness_verified": True,
                      "fixed_order64_exact_defect": str(Fraction(numerator, 64 ** 2))}, indent=2))
