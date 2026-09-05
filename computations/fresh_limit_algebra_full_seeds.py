"""Numerical PSD-majorant screening of full diagonal completions.

These are numerical SDP optima, not exact upper certificates.
"""
import itertools
import json
import numpy as np
import cvxpy as cp


def screen(A):
    n = len(A)
    X = np.array(list(itertools.product((-1, 1), repeat=n)))
    D = cp.Variable((n, n), symmetric=True)
    t = cp.Variable()
    Ap = cp.Parameter((n, n), symmetric=True)
    constraints = [D - Ap >> 0, D + Ap >> 0]
    constraints.extend(cp.sum(cp.multiply(np.outer(x, x), D)) <= 2*t for x in X)
    problem = cp.Problem(cp.Minimize(t), constraints)
    records = []
    for diagonal in itertools.product((-1, 1), repeat=n):
        B = A + np.diag(diagonal)
        Ap.value = B
        value = problem.solve(solver='CLARABEL', tol_gap_abs=1e-9, tol_feas=1e-9,
                              tol_gap_rel=1e-9, max_iter=200)
        q = float(np.max(np.abs(np.sum((X @ B) * X, axis=1))) / 2)
        records.append({'diagonal': diagonal, 'T_numeric': value, 'Q_exact': q})
    records.sort(key=lambda x: x['T_numeric'])
    return {'n': n, 'half_floor': n**1.5/2, 'best': records[:8], 'worst': records[-1]}


A5 = np.array([[0,-1,1,-1,1],[-1,0,-1,1,1],[1,-1,0,1,-1],
               [-1,1,1,0,-1],[1,1,-1,-1,0]])
A6 = np.ones((6, 6), dtype=int) - np.eye(6, dtype=int)
for bit, (i,j) in enumerate(itertools.combinations(range(1, 6), 2)):
    A6[i,j] = A6[j,i] = 1 - 2*((220 >> bit) & 1)
print(json.dumps(screen(A5)))
print(json.dumps(screen(A6)))
