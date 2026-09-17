"""Exact CP-SAT search for a three-cell reflection quotient of H144.

stdout only; a found coloring is verified with integer matrix products.
"""

import json
import numpy as np
from ortools.sat.python import cp_model

p = 11
sq = {i*i % p for i in range(1,p)}
F = np.ones((12,12),dtype=np.int64)
F[1:,0] = -1
for i in range(p):
    for j in range(p):
        F[i+1,j+1] = 1 if i == j or (j-i)%p in sq else -1
assert np.array_equal(F.T @ F,12*np.eye(12,dtype=np.int64))

model = cp_model.CpModel()
x = [[[model.NewBoolVar(f'x_{i}_{j}_{c}') for c in range(3)]
      for j in range(12)] for i in range(12)]
for i in range(12):
    for j in range(12):
        model.AddExactlyOne(x[i][j])
for j in range(12):
    for c in range(3):
        model.Add(sum(x[i][j][c] for i in range(12)) == 4)
for j in range(12):
    for k in range(j+1,12):
        for c in range(3):
            model.Add(sum(int(F[i,j]*F[i,k])*(x[i][j][c]+x[i][k][c])
                          for i in range(12)) == 0)
model.Add(x[0][0][0] == 1)
solver = cp_model.CpSolver()
solver.parameters.max_time_in_seconds = 240
solver.parameters.num_search_workers = 4
solver.parameters.random_seed = 20260906
status = solver.Solve(model)
print(solver.StatusName(status), flush=True)
if status in (cp_model.FEASIBLE,cp_model.OPTIMAL):
    colors = np.array([[next(c for c in range(3) if solver.Value(x[i][j][c]))
                        for j in range(12)] for i in range(12)])
    for c in range(3):
        Z = F*(colors == c)
        assert np.array_equal(Z.T@F+F.T@Z,8*np.eye(12,dtype=np.int64))
    print(json.dumps(colors.tolist()), flush=True)
print(solver.ResponseStats(), flush=True)
