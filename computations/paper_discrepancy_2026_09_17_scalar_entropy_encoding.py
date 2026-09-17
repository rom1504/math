"""Exact importance-encoding identities for isotropic Walsh sign laws."""

from __future__ import annotations

from fractions import Fraction
import json
import math

import numpy as np

from paper_discrepancy_2026_09_17_scalar_entropy_counterexample import graph_words, parity


def check_word(walsh, word, budget):
    n = len(word)
    projections = walsh @ word
    response = Fraction(int(np.abs(projections).sum()), n)
    assert response <= budget
    reference_mass = Fraction(1, 2*n)
    signed_atoms, density = [], []
    for atom, projection in zip(walsh, projections):
        for sign in (-1, 1):
            signed_atoms.append(sign*atom)
            density.append(1-response/budget
                           +Fraction(2*max(sign*int(projection), 0), budget))
    assert min(density) >= 0
    assert sum(density)*reference_mass == 1
    for coordinate in range(n):
        mean = sum(q*int(atom[coordinate]) for q, atom in zip(density, signed_atoms))*reference_mass
        assert mean == Fraction(int(word[coordinate]), budget)
    second = sum(q*q for q in density)*reference_mass
    assert second == 1+(2*n-response**2)/budget**2
    divergence = sum(float(q)*math.log(float(q)) for q in density if q)*float(reference_mass)
    assert divergence <= math.log(float(second))+1e-12
    assert math.log(float(second)) <= math.log(1+2*n/budget**2)+1e-12
    return {"response": str(response), "density_second_moment": str(second),
            "KL": divergence, "KL_upper": math.log(1+2*n/budget**2)}


def main():
    reports = []
    for r in range(1, 4):
        n = 4**r
        walsh = np.asarray([[(-1)**parity(i & j) for j in range(n)]
                            for i in range(n)], dtype=np.int64)
        assert np.array_equal(walsh.T @ walsh, n*np.eye(n, dtype=np.int64))
        count = 0
        for word in graph_words(r):
            result = check_word(walsh, word, 3)
            count += 1
        reports.append({"r": r, "dimension": n, "checked_words": count,
                        "last_word": result})
    print(json.dumps({"status": "all importance-encoding checks passed", "cases": reports}, indent=2))


if __name__ == "__main__":
    main()
