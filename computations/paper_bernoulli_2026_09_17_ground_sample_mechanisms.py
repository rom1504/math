"""Exact physical ground-sample products and majority; no LP solver.

Independent global centering of each sample is retained in majority.
All response numerators, covariances and support masses are integers.
"""

from fractions import Fraction as F
import importlib.util
from itertools import combinations, product
import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def fwht(values):
    out = np.asarray(values, dtype=np.int64).copy()
    step = 1
    while step < len(out):
        blocks = out.reshape(-1, 2 * step)
        left = blocks[:, :step].copy()
        right = blocks[:, step:].copy()
        blocks[:, :step] = left + right
        blocks[:, step:] = left - right
        step *= 2
    return out


def inverse(values):
    out = fwht(values)
    assert np.all(out % len(out) == 0)
    return out // len(out)


def mask_of(word):
    return sum((int(v) < 0) << i for i, v in enumerate(word[1:]))


def main():
    spec = importlib.util.spec_from_file_location(
        "stored", ROOT / "computations/transfer_adversary_minimizer_isotropy_2026_09_06.py")
    stored = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(stored)
    reports = []
    for name, matrix, provenance in stored.cases():
        n = len(matrix)
        if n not in (6, 12, 14):
            continue
        a = np.asarray(matrix, dtype=np.int64)
        words = np.asarray([(1,) + s for s in product((-1, 1), repeat=n-1)], dtype=np.int64)
        edges = list(combinations(range(n), 2))
        energies = sum(a[i, j] * words[:, i] * words[:, j] for i, j in edges)
        cap = int(max(abs(energies)))
        assert cap == provenance["minimum_cap_imported"]
        code = words[abs(energies) == cap]
        masks = np.asarray([mask_of(x) for x in code], dtype=np.int64)
        length = 1 << (n-1)
        initial = np.bincount(masks, minlength=length)
        transformed = fwht(initial)
        kernel = np.asarray([abs(n - 2 * bin(i).count("1")) for i in range(length)], dtype=np.int64)
        kernel_transform = fwht(kernel)
        distributions = [("one_ground_sample", initial)]
        for power in (2, 3, 4):
            distributions.append((f"product_{power}", inverse(transformed**power)))
        # Fix the first sample's global sign by overall-output symmetry,
        # but center the SECOND and THIRD samples independently.
        full_mask = (1 << n) - 1
        projective_mask = length - 1
        first_bit = length
        centered_samples = np.concatenate((masks, masks ^ full_mask))
        b = centered_samples[:, None]
        c = centered_samples[None, :]
        majority_counts = np.zeros(length, dtype=np.int64)
        for x in masks:
            majority = (x & b) | (x & c) | (b & c)
            normalized = np.where((majority & first_bit) != 0,
                                  majority ^ full_mask, majority) & projective_mask
            majority_counts += np.bincount(normalized.ravel(), minlength=length)
        distributions.append(("centered_majority_3", majority_counts))
        for mechanism, counts in distributions:
            assert min(counts) >= 0
            total = int(sum(counts))
            transform = fwht(counts)
            for i in range(n-1):
                assert transform[1 << i] == 0
                for j in range(i+1, n-1):
                    assert transform[(1 << i) | (1 << j)] == 0
            response = inverse(transform * kernel_transform)
            ground_response = response[masks]
            optimal = {6: F(5, 3), 12: F(9, 5), 14: F(98, 39)}[n]
            minimum = F(int(min(ground_response)), total)
            maximum = F(int(max(ground_response)), total)
            average = F(int(sum(ground_response)), total * len(masks))
            assert average >= optimal
            reports.append({
                "case": name, "n": n, "cap": cap,
                "projective_ground_size": len(code), "mechanism": mechanism,
                "sample_count_denominator": total,
                "physical_support_size": int(np.count_nonzero(counts)),
                "exact_isotropy": True,
                "exact_ground_support_mass": str(F(int(sum(counts[masks])), total)),
                "exact_min_ground_response": str(minimum),
                "exact_max_ground_response": str(maximum),
                "exact_average_ground_response": str(average),
                "exact_optimum": str(optimal),
                "exact_max_over_optimum": str(maximum / optimal),
                "normalized_max_response": float(maximum) / n**.5,
                "normalized_target": 3 * cap / (2 * n**1.5),
            })
    result = {"status": "PASS exact physical-law covariance and response",
              "cases": reports,
              "scope": "Explicit finite mechanism tests, not asymptotic optimizer conclusions."}
    output = ROOT / "tmp/paper_portfolio_2026_09_17/bernoulli/ground_sample_mechanisms.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
