"""Independent, solver-free replay of full physical-column game witnesses.

Pass --witness to either the live JSON or its preserved archive copy.
No numerical optimizer, response-game implementation, or floating tolerances
are used in the mathematical checks. Stored matrix optimality labels remain
imported; their actual caps are exhaustively recomputed.
"""

import argparse
from fractions import Fraction
import importlib.util
import itertools
import json
import math
from pathlib import Path

import numpy as np


def common_integers(values):
    denominator = 1
    for value in values:
        denominator = math.lcm(denominator, value.denominator)
    return denominator, np.array(
        [value.numerator * (denominator // value.denominator) for value in values],
        dtype=object,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--witness", required=True, type=Path)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    spec = importlib.util.spec_from_file_location(
        "matrix_inputs_only", root / "computations/transfer_adversary_minimizer_isotropy_2026_09_06.py"
    )
    inputs = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(inputs)
    matrices = {name: np.asarray(matrix, dtype=np.int64)
                for name, matrix, _ in inputs.cases()}
    witness = json.loads(args.witness.read_text())
    reports = []
    for record in witness["cases"]:
        a = matrices[record["case"]]
        n = len(a)
        assert n == record["n"]
        assert np.array_equal(a, a.T) and not np.diag(a).any()
        edges = list(itertools.combinations(range(n), 2))
        assert all(abs(int(a[i, j])) == 1 for i, j in edges)
        spins = np.array([(1,) + tail for tail in itertools.product((-1, 1), repeat=n-1)],
                         dtype=np.int64)
        features = np.array([spins[:, i] * spins[:, j] for i, j in edges], dtype=np.int64).T
        energy = features @ np.array([a[i, j] for i, j in edges], dtype=np.int64)
        cap = int(np.abs(energy).max())
        assert cap == record["cap"]
        code = spins[cap - np.abs(energy) <= record["deficit_window"]]
        assert len(code) == record["projective_code_size"]
        code_set = {tuple(x) for x in code}

        primal = record["primal_law"]
        weights = [Fraction(item["weight"]) for item in primal]
        assert all(w >= 0 for w in weights) and sum(weights) == 1
        words = np.array([item["word"] for item in primal], dtype=np.int64)
        assert words.shape == (len(primal), n) and np.all(np.abs(words) == 1)
        pd, pn = common_integers(weights)
        responses = np.abs(code @ words.T).astype(object) @ pn
        upper = Fraction(int(max(responses)), pd)
        assert upper == Fraction(record["upper_exact"])
        if record["isotropy_required"]:
            correlations = np.array([words[:, i] * words[:, j] for i, j in edges], dtype=object) @ pn
            assert all(value == 0 for value in correlations)

        dual = record["dual_code_law"]
        pi = [Fraction(item["weight"]) for item in dual]
        assert all(value >= 0 for value in pi) and sum(pi) <= 1
        dual_words = np.array([item["word"] for item in dual], dtype=np.int64).reshape((-1, n))
        assert all(tuple(x) in code_set for x in dual_words)
        coefficients = [Fraction(s) for s in record["dual_equality_coefficients"]]
        expected = 1 + (len(edges) if record["isotropy_required"] else 0)
        assert len(coefficients) == expected
        dd, dn = common_integers(coefficients + pi)
        lhs = np.abs(spins @ dual_words.T).astype(object) @ dn[len(coefficients):]
        rhs = np.full(len(spins), dn[0], dtype=object)
        if record["isotropy_required"]:
            rhs += features.astype(object) @ dn[1:len(coefficients)]
        margin = min(lhs - rhs)
        assert margin >= 0
        lower = coefficients[0]
        assert lower == Fraction(record["lower_exact"]) and lower <= upper
        reports.append({"case": record["case"], "window": record["deficit_window"],
                        "isotropic": record["isotropy_required"], "lower": str(lower),
                        "upper": str(upper), "all_columns_checked": len(spins),
                        "all_code_words_checked": len(code),
                        "minimum_dual_margin": str(Fraction(int(margin), dd))})
        print(json.dumps({"passed": len(reports), "case": record["case"],
                          "window": record["deficit_window"],
                          "isotropic": record["isotropy_required"]}), flush=True)
    print(json.dumps({"status": "PASS solver-free exact rational replay", "games": reports,
                      "global_signing_minimality_not_reproved": True}, indent=2))


if __name__ == "__main__":
    main()
