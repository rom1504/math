"""Exact all-signed-matching dual test on finite actual ground codes.

Law: one global fair pole, independent fair pair signs, prescribed
signed perfect matching. All such laws are physical and isotropic.
"""

from collections import Counter
from fractions import Fraction as F
import importlib.util
from itertools import combinations, product
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def matchings(indices):
    if not indices:
        yield ()
        return
    first = indices[0]
    for j in range(1, len(indices)):
        other = indices[j]
        remaining = indices[1:j] + indices[j+1:]
        for tail in matchings(remaining):
            yield ((first, other),) + tail


def main():
    spec = importlib.util.spec_from_file_location(
        "stored", ROOT / "computations/transfer_adversary_minimizer_isotropy_2026_09_06.py")
    stored = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(stored)
    reports = []
    for name, matrix, provenance in stored.cases():
        n = len(matrix)
        if n not in (6, 8, 12):
            continue
        m = n // 2
        a = np.asarray(matrix, dtype=np.int64)
        words = np.asarray([(1,) + s for s in product((-1, 1), repeat=n-1)], dtype=np.int64)
        energies = sum(a[i, j] * words[:, i] * words[:, j]
                       for i, j in combinations(range(n), 2))
        cap = int(max(abs(energies)))
        assert cap == provenance["minimum_cap_imported"]
        code = words[abs(energies) == cap]
        absolute_means = [F(sum(math.comb(k, j) * abs(k-2*j) for j in range(k+1)), 2**k)
                          for k in range(m+1)]
        costs = [absolute_means[k] + absolute_means[m-k] for k in range(m+1)]
        den = math.lcm(*(c.denominator for c in costs))
        lookup = np.asarray([int(c * den) for c in costs], dtype=np.int64)
        # Simultaneous reversal of all edge signs swaps the two poles.
        signs = np.asarray([(1,) + s for s in product((-1, 1), repeat=m-1)], dtype=np.int64)
        best_total = 10**18
        best_maximum = 10**18
        best_witness = None
        best_histogram = None
        best_response_sum = np.zeros(len(code), dtype=np.int64)
        count = matching_count = minimizers = 0
        for pairs in matchings(tuple(range(n))):
            edge_features = np.asarray([code[:, i] * code[:, j] for i, j in pairs], dtype=np.int64).T
            counts = (m + edge_features @ signs.T) // 2
            response = lookup[counts]
            totals = response.sum(axis=0)
            maxima = response.max(axis=0)
            candidate = int(totals.min())
            if candidate < best_total:
                best_total = candidate
                minimizers = int(np.count_nonzero(totals == candidate))
                best_response_sum = response[:, totals == candidate].sum(axis=1)
                index = int(np.argmin(totals))
                best_witness = {"pairs": [list(e) for e in pairs],
                                "edge_signs": signs[index].tolist()}
                best_histogram = {str(F(int(v), den)): int(k)
                                  for v, k in Counter(map(int, response[:, index])).items()}
            elif candidate == best_total:
                minimizers += int(np.count_nonzero(totals == candidate))
                best_response_sum += response[:, totals == candidate].sum(axis=1)
            best_maximum = min(best_maximum, int(maxima.min()))
            count += len(signs)
            matching_count += 1
        expected_matchings = math.prod(range(1, n, 2))
        assert matching_count == expected_matchings
        lower = F(best_total, den * len(code))
        uniform_best_min = F(int(min(best_response_sum)), den * minimizers)
        uniform_best_max = F(int(max(best_response_sum)), den * minimizers)
        single_upper = F(best_maximum, den)
        target = F(3 * cap, 2 * n)
        reports.append({
            "case": name, "n": n, "cap": cap, "ground_code_size": len(code),
            "perfect_matchings_exhausted": matching_count,
            "signed_matching_laws_exhausted_up_to_pole_exchange": count,
            "exact_uniform_ground_dual_lower_for_ALL_matching_mixtures": str(lower),
            "exact_best_single_matching_max_response": str(single_upper),
            "exact_target_3Q_over_2n": str(target),
            "exact_dual_minus_target": str(lower-target),
            "uniform_dual_optimal_law_count": minimizers,
            "exact_uniform_mixture_of_dual_minimizers_min_response": str(uniform_best_min),
            "exact_uniform_mixture_of_dual_minimizers_max_response": str(uniform_best_max),
            "exact_matching_class_game_certified": str(lower) if uniform_best_max == lower else None,
            "one_minimizing_matching": best_witness,
            "its_exact_ground_response_histogram": best_histogram,
        })
    result = {"status": "PASS exhaustive integer signed-matching response test",
              "global_minimality_imported": True, "cases": reports,
              "scope": "Finite exact class obstruction; no asymptotic impossibility for arbitrary physical laws."}
    output = ROOT / "tmp/paper_portfolio_2026_09_17/bernoulli/signed_matching_ground_game.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
