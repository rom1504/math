#!/usr/bin/env python3
"""Exact partial-Hadamard saturation checks; no Boolean-energy assertion."""

import argparse
import json
from fractions import Fraction
from pathlib import Path
from random import Random


def hadamard(order):
    return [[1 if bin(i & j).count("1") % 2 == 0 else -1
             for j in range(order)] for i in range(order)]


def nullspace(matrix, columns):
    a = [[Fraction(v) for v in row] for row in matrix]
    pivots = []
    row = 0
    for col in range(columns):
        pivot = next((r for r in range(row, len(a)) if a[r][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        scale = a[row][col]
        a[row] = [v / scale for v in a[row]]
        for r in range(len(a)):
            if r != row and a[r][col]:
                scale = a[r][col]
                a[r] = [v - scale * w for v, w in zip(a[r], a[row])]
        pivots.append(col)
        row += 1
        if row == len(a):
            break
    basis = []
    for free in set(range(columns)) - set(pivots):
        vector = [Fraction(0) for _ in range(columns)]
        vector[free] = Fraction(1)
        for r, pivot in enumerate(pivots):
            vector[pivot] = -a[r][free]
        basis.append(vector)
    return basis


def check_frame(order, selected, core):
    full = hadamard(order)
    f = [full[i] for i in selected]
    q, s = len(selected), len(core)
    assert 0 not in core and s > order - q
    assert all(sum(f[i][k] * f[j][k] for k in range(order))
               == (order if i == j else 0)
               for i in range(q) for j in range(q))
    complement = sorted(set(range(order)) - set(core))
    basis = nullspace([[f[i][j] for i in range(q)] for j in complement], q)
    assert len(basis) >= q + s - order > 0
    # B=q*P*F_C is integer, so its Gram has eigenvalue q^2*order.
    b = [[q * f[i][j] - sum(f[r][j] for r in range(q)) for j in core]
         for i in range(q)]
    assert all(sum(b[i][j] for i in range(q)) == 0 for j in range(s))
    for vector in basis:
        assert sum(vector) == 0
        response = [sum(b[i][j] * vector[i] for i in range(q)) for j in range(s)]
        assert all(sum(b[i][j] * response[j] for j in range(s))
                   == q * q * order * vector[i] for i in range(q))
    return {"m": order, "q": q, "s": s, "selected_rows": selected,
            "active_core": core, "forced_multiplicity": q + s - order,
            "exact_complement_kernel_dimension": len(basis)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    rng = Random(2026091645)
    checks = []
    for order in (4, 8, 16, 32):
        for q in sorted({order // 2, 3 * order // 4, order - 2, order}):
            for _ in range(3):
                selected = sorted(rng.sample(range(order), q))
                for s in sorted({order - q + 1, order - 1}):
                    if s <= order - 1:
                        core = sorted(rng.sample(range(1, order), s))
                        checks.append(check_frame(order, selected, core))
    # A boundary case has no forced saturation: m=4,q=3,s=1.
    # The chosen non-DC column is (1,-1,1), whose centered norm²=8/3<4.
    boundary_column = [hadamard(4)[i][1] for i in (0, 1, 2)]
    assert sum(v * v for v in boundary_column) - Fraction(sum(boundary_column)**2, 3) == Fraction(8, 3)
    report = {"status": "PASS", "exact_cases": len(checks), "checks": checks,
              "boundary_example": {"m": 4, "q": 3, "s": 1,
                                   "centered_squared_norm": "8/3 < 4"},
              "scope": "Exact rational frame identities, not Boolean energy witnesses; repair transfer is analytic."}
    rendered = json.dumps(report, indent=2, sort_keys=True)
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered + "\n")
    print(rendered)


if __name__ == "__main__":
    main()
