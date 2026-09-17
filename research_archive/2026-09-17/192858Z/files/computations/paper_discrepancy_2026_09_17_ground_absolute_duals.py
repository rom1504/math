"""Finite unrestricted absolute-response dual certificates for stored minimizers."""

from __future__ import annotations

from fractions import Fraction
import argparse
import importlib.util
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog


ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--full", action="store_true", help="Print all rational certificate atoms.")
    arguments = parser.parse_args()
    source = ROOT / "computations/transfer_adversary_minimizer_isotropy_2026_09_06.py"
    spec = importlib.util.spec_from_file_location("stored_minimizers", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    reports = []
    for name, matrix, provenance in module.cases():
        n = len(matrix)
        if n % 2 or n > 14:
            continue
        a = np.asarray(matrix, dtype=np.int64)
        words = np.asarray([(1,)+row for row in itertools.product((-1, 1), repeat=n-1)],
                           dtype=np.int64)
        energies = np.einsum("bi,ij,bj->b", words, a, words)//2
        cap = int(np.abs(energies).max())
        assert cap == provenance["minimum_cap_imported"]
        ground = words[np.abs(energies) == cap]
        overlaps = np.abs(words @ ground.T)
        uniform_min = Fraction(int(overlaps.sum(axis=1).min()), len(ground))
        count = len(ground)
        objective = np.zeros(count+1)
        objective[-1] = -1
        upper = np.column_stack((-overlaps, np.ones(len(words))))
        equality = np.zeros((1,count+1))
        equality[0,:count] = 1
        result = linprog(objective, A_ub=upper, b_ub=np.zeros(len(words)),
                         A_eq=equality, b_eq=[1], bounds=[(0,None)]*(count+1),
                         method="highs")
        assert result.success, result.message
        denominator = 10000
        target = denominator*result.x[:count]
        numerator = np.floor(np.maximum(target,0)).astype(np.int64)
        remainder = denominator-int(numerator.sum())
        order = np.argsort(-(target-numerator))
        numerator[order[:remainder]] += 1
        assert int(numerator.sum()) == denominator and np.all(numerator >= 0)
        integer_responses = overlaps @ numerator
        certified = Fraction(int(integer_responses.min()),denominator)
        support = np.flatnonzero(numerator)
        report = {"case":name,"n":n,"cap":cap,
                        "ground_projective_words":count,
                        "uniform_ground_dual":str(uniform_min),
                        "unrestricted_game_optimum_numeric":float(result.x[-1]),
                        "rational_certificate":str(certified),
                        "normalized_certificate":float(certified)/np.sqrt(n),
                        "support_size":len(support),
                        "certificate_denominator":denominator,
                        "support_words":[ground[i].tolist() for i in support],
                        "support_numerators":numerator[support].tolist()}
        if not arguments.full:
            del report["support_words"]
            del report["support_numerators"]
        reports.append(report)
    print(json.dumps({"status":"all rational all-sign-query certificates passed",
                      "scope":"Global minimizer labels imported; dual lower bounds exhaust all sign queries.",
                      "cases":reports},indent=2))


if __name__ == "__main__":
    main()
