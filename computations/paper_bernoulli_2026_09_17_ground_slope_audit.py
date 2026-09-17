"""Solver-free replay of physical ground response versus 3Q/(2n).

Caps and full primal/dual inequalities are exact. The stored matrices'
global-minimizer labels remain imported provenance, not newly proved.
"""

from fractions import Fraction as F
import importlib.util
from itertools import combinations, product
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]


def integers(values):
    den = math.lcm(*(v.denominator for v in values))
    return den, np.asarray([v.numerator * (den // v.denominator)
                            for v in values], dtype=object)


def main():
    spec = importlib.util.spec_from_file_location(
        "stored_minimizers", ROOT / "computations/transfer_adversary_minimizer_isotropy_2026_09_06.py")
    stored = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(stored)
    inputs = json.loads((ROOT / "tmp/paper_portfolio_2026_09_17/discrepancy/full_column_game_certificates.json").read_text())
    matrices = {name: (matrix, provenance)
                for name, matrix, provenance in stored.cases()}
    reports = []
    checked_columns = checked_queries = 0
    for certificate in inputs["cases"]:
        n = certificate["n"]
        window = certificate["deficit_window"]
        if window != 0 and not (n in (12, 14) and window == 2):
            continue
        name = certificate["case"]
        matrix, provenance = matrices[name]
        a = np.asarray(matrix, dtype=np.int64)
        words = np.asarray([(1,) + tail for tail in product((-1, 1), repeat=n-1)], dtype=np.int64)
        edges = list(combinations(range(n), 2))
        features = np.asarray([words[:, i] * words[:, j] for i, j in edges], dtype=np.int64).T
        energies = features @ np.asarray([a[i, j] for i, j in edges], dtype=np.int64)
        cap = int(max(abs(energies)))
        assert cap == certificate["cap"] == provenance["minimum_cap_imported"]
        code = words[cap - abs(energies) <= window]
        assert len(code) == certificate["projective_code_size"]
        primal_words = np.asarray([row["word"] for row in certificate["primal_law"]], dtype=np.int64)
        primal_weights = [F(row["weight"]) for row in certificate["primal_law"]]
        assert all(w >= 0 for w in primal_weights) and sum(primal_weights) == 1
        den, numerator = integers(primal_weights)
        responses = abs(code @ primal_words.T).astype(object) @ numerator
        upper = F(int(max(responses)), den)
        assert upper == F(certificate["upper_exact"])
        if certificate["isotropy_required"]:
            covariance = (primal_words.T.astype(object) * numerator) @ primal_words
            assert np.array_equal(covariance, den * np.eye(n, dtype=object))
        dual_rows = certificate["dual_code_law"]
        dual_words = np.asarray([row["word"] for row in dual_rows], dtype=np.int64)
        dual_weights = [F(row["weight"]) for row in dual_rows]
        assert all(w >= 0 for w in dual_weights) and sum(dual_weights) <= 1
        code_set = set(map(tuple, code))
        assert all(tuple(word) in code_set for word in dual_words)
        coefficients = [F(v) for v in certificate["dual_equality_coefficients"]]
        den, numerator = integers(coefficients + dual_weights)
        if dual_rows:
            lhs = abs(words @ dual_words.T).astype(object) @ numerator[len(coefficients):]
        else:
            lhs = np.zeros(len(words), dtype=object)
        rhs = np.full(len(words), numerator[0], dtype=object)
        if len(coefficients) > 1:
            assert len(coefficients) == len(edges) + 1
            rhs += features.astype(object) @ numerator[1:len(coefficients)]
        assert min(lhs - rhs) >= 0
        lower = coefficients[0]
        assert lower == upper == F(certificate["lower_exact"])
        target = F(3 * cap, 2 * n)
        reports.append({
            "case": name, "n": n, "cap": cap, "window": window,
            "isotropy_required": certificate["isotropy_required"],
            "exact_game": str(upper), "exact_target_3Q_over_2n": str(target),
            "exact_game_over_target": str(upper / target),
            "normalized_cap": cap / n**1.5,
            "normalized_response": float(upper) / math.sqrt(n),
            "normalized_target": float(target) / math.sqrt(n),
            "literal_target_holds": upper <= target,
        })
        checked_columns += len(words)
        checked_queries += len(code)
    result = {
        "status": "PASS exact all-column primal/dual replay",
        "global_minimality_imported": True,
        "all_column_lower_checks": checked_columns,
        "all_query_upper_checks": checked_queries,
        "cases": reports,
        "scope": "Finite target test only; no asymptotic optimizer obstruction follows.",
    }
    output = ROOT / "tmp/paper_portfolio_2026_09_17/bernoulli/ground_slope_audit.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
