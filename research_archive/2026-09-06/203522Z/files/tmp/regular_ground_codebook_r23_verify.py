#!/usr/bin/env python3
"""Exact finite row-field audit for the regular-ground codebooks."""

from fractions import Fraction
import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "selector_codebook_r22", ROOT / "tmp" / "selector_codebook_r22.py"
)
SC = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(SC)


def audit(name, A, mass, _expected_counts, _expected_z):
    patterns = SC.coverage_patterns(A)
    parent_q, _ = SC.PC.matrix_norm(A)
    profiles = []
    for pattern in mass:
        representatives = [
            (t, x)
            for t, x, parent_score in patterns[pattern]
            if parent_score == parent_q
        ]
        assert representatives
        t, x = representatives[0]
        row = [
            sum(t * A[i][j] * x[i] * x[j] for j in range(len(A)) if j != i)
            for i in range(len(A))
        ]
        assert min(row) >= 0
        assert sum(row) == 2 * parent_q
        assert sum(value * value for value in row) <= max(row) * 2 * parent_q
        profiles.append(tuple(row))
    print(
        name,
        "support_size=", len(mass),
        "ordered_Q=", 2 * parent_q,
        "max_row=", max(max(row) for row in profiles),
        "max_R2=", max(sum(value * value for value in row) for row in profiles),
    )


def main():
    for name, data in SC.CASES.items():
        audit(name, *data)
    print("all regular-ground finite audits passed")


if __name__ == "__main__":
    main()
