#!/usr/bin/env python3
"""Exact small certificate disproving recovery at rotation velocity."""
import sympy as sp


def main():
    values = []
    for mask in range(8):
        x, y, z = [1 - 2 * ((mask >> j) & 1) for j in range(3)]
        value = (1 - x*z + sp.I*y*(x+z)) / 2
        assert sp.expand(value * sp.conjugate(value)) == 1
        values.append(value)
        if x == z:
            assert value == sp.I*x*y
    J = sp.zeros(8)
    for x in range(8):
        for j in range(3):
            J[x, x ^ (1 << j)] = sp.Rational(1, 2)
    U = sp.diag(*values)
    A = (U.conjugate().T * J * U - J).applyfunc(sp.expand)
    assert A == A.conjugate().T
    assert (A**3 - 3*A).applyfunc(sp.expand) == sp.zeros(8)
    assert sp.trace(A**2) == 12
    polynomial = sp.factor(A.charpoly().as_expr())
    print("Values:", values)
    print("Characteristic polynomial:", polynomial)
    print("PASS: A^3=3A and tr(A^2)=12, so v(f)=sqrt(3).")
    print("PASS: x3=x1 restricts f to i*x1*x2; tensor uniform degree rate=2.")


if __name__ == "__main__":
    main()
