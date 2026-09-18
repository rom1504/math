#!/usr/bin/env python3
"""Exact characteristic polynomial certificate for the five-bit velocity.

The bipartite block C of 2*(f^* J f-J) has Eisenstein integer entries.
Eigenvalues of C C^* are four times the squared positive velocities.
All polynomial and root-isolation operations below are exact.
"""
import importlib.util
from pathlib import Path

import sympy as sp
from sympy.polys.matrices import DomainMatrix


def main():
    path = Path(__file__).with_name("bh_boundedness_barrier_2026_09_18_tangent.py")
    spec = importlib.util.spec_from_file_location("exact_tangent", path)
    t = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(t)
    roots = [(1, 0), (0, 1), (-1, 1), (-1, 0), (0, -1), (1, -1)]
    masks = [0, 1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 16, 17, 18, 20, 24]
    powers = [0, 2, 0, 5, 1, 4, 2, 4, 3, 1, 4, 1, 2, 2, 3, 2]
    f4 = dict(zip(masks, [roots[j] for j in powers]))
    fv = []
    for x in range(32):
        u, v = t.evaluate(f4, x)
        assert u % 4 == v % 4 == 0
        fv.append((u // 4, v // 4))
        assert fv[-1] in roots
    evens = [x for x in range(32) if bin(x).count("1") % 2 == 0]
    odds = [x for x in range(32) if bin(x).count("1") % 2 == 1]
    C = []
    for x in evens:
        row = []
        for y in odds:
            value = (0, 0)
            if bin(x ^ y).count("1") == 1:
                value = t.add(t.mul(t.conj(fv[x]), fv[y]), (-1, 0))
            row.append(value)
        C.append(row)
    D = []
    zeta = (1 + sp.sqrt(-3)) / 2
    for i in range(16):
        row = []
        for j in range(16):
            value = (0, 0)
            for k in range(16):
                value = t.add(value, t.mul(C[i][k], t.conj(C[j][k])))
            row.append(value[0] + value[1] * zeta)
        D.append(row)
    field = sp.QQ.algebraic_field(sp.sqrt(-3))
    mat = DomainMatrix.from_list_sympy(16, 16, D).convert_to(field)
    coefficients = [field.to_sympy(c).expand() for c in mat.charpoly()]
    assert all(c.is_Rational for c in coefficients)
    x = sp.Symbol("x")
    polynomial = sp.Poly.from_list(coefficients, x, domain=sp.QQ)
    print("Characteristic polynomial of C C*:", polynomial.as_expr(), flush=True)
    print("Factorization:", sp.factor_list(polynomial.as_expr()), flush=True)
    intervals = polynomial.intervals(eps=sp.Rational(1, 10**10))
    print("Exact isolating intervals and multiplicities:", intervals, flush=True)
    last_interval, multiplicity = intervals[-1]
    assert multiplicity == 1
    assert all(interval[1] < last_interval[0] for interval, _ in intervals[:-1])
    assert last_interval[0] > sp.Rational(1392, 100)
    assert last_interval[1] < sp.Rational(1393, 100)
    print("PASS: largest positive velocity is simple; 1.865 < v(f) < 1.867.")


if __name__ == "__main__":
    main()
