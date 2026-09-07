#!/usr/bin/env python3
"""Evaluate (10.426) on the strict nine-type weighted obstruction."""

from fractions import Fraction as F
import importlib.util
from pathlib import Path
import sys


HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "minimax_allocation_r8", HERE / "minimax_allocation_r8.py"
)
mm = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mm
spec.loader.exec_module(mm)


S = set(range(4))
B = {4, 5, 6}


def matrix_capacity():
    a = [[F(0) for _ in range(9)] for _ in range(9)]
    for i in range(9):
        for j in range(i + 1, 9):
            if i in S and j in S:
                z = 54
            elif i in S and j in B:
                z = 28
            elif i in S and j == 7:
                z = 55
            elif i in S and j == 8:
                z = 29
            elif i in B and j in B:
                z = 83
            elif i in B and j in {7, 8}:
                z = -55
            elif (i, j) == (7, 8):
                z = 113
            else:
                raise AssertionError((i, j))
            a[i][j] = a[j][i] = F(z, 672)
    return tuple(tuple(row) for row in a)


def matrix_strictified():
    u, v, c, d = set(range(3)), set(range(3, 5)), set(range(5, 8)), {8}
    a = [[F(0) for _ in range(9)] for _ in range(9)]
    for i in range(9):
        for j in range(i + 1, 9):
            if i in u and j in u:
                z = F(49, 400)
            elif i in u and j in v:
                z = F(-49, 600)
            elif i in v and j in v:
                z = F(49, 300)
            elif i in c and j in c:
                z = F(49, 400)
            elif i in c and j in d:
                z = F(49, 1200)
            elif i in u and j in c:
                z = F(251, 6000)
            elif i in u and j in d:
                z = F(251, 6000)
            elif i in v and j in c:
                z = F(31, 375)
            elif i in v and j in d:
                z = F(1, 1000)
            else:
                raise AssertionError((i, j))
            a[i][j] = a[j][i] = z
    return tuple(tuple(row) for row in a)


def evaluate(label, a):
    print("LABEL", label)
    P, N, ps, ns = mm.extrema(a)
    print("P,N,R,endpoints", P, N, P + N, len(ps), len(ns))
    sigs = mm.tree_signatures(a)
    print("signatures", len(sigs))
    rows = []
    for k, sig in enumerate(sigs):
        root = mm.node_from_signature(sig)
        ans, _, _ = mm.primal_lp(root, P + N)
        rows.append((F(float(ans.fun)).limit_denominator(10**8), k, sig))
    rows.sort(key=lambda z: (z[0], repr(z[2])))
    print("K min/max", rows[0][0], rows[-1][0])
    print("min index/signature", rows[0][1], rows[0][2])
    root = mm.node_from_signature(rows[0][2])
    print("exact primal", mm.exact_primal_certificate(root, P + N)[0])
    print("exact dual", mm.exact_dual_certificate(root, P + N)[0])
    mm.summarize(root, P + N)


def main():
    evaluate("capacity", matrix_capacity())
    evaluate("strictified", matrix_strictified())


if __name__ == "__main__":
    main()
