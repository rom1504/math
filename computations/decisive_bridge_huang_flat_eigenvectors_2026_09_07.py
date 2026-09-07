"""Exact integer Boolean eigenvector search and multiplicative verification."""
import argparse
import json
import math
import numpy as np
from ortools.sat.python import cp_model


def apply_h(v, d):
    if d == 0:
        return np.zeros_like(v)
    half = len(v)//2
    return np.r_[apply_h(v[:half], d-1)+v[half:],
                 v[:half]-apply_h(v[half:], d-1)]


def neighbors(d, vertex):
    for j in range(d):
        # Recurrence orders the highest bit first; prior higher bits give sign.
        sign = (-1)**(bin(vertex >> (j+1)).count("1"))
        yield vertex ^ (1 << j), sign


def solve(d, seconds):
    root = math.isqrt(d)
    assert root*root == d
    n=1<<d
    model=cp_model.CpModel()
    b=[model.new_bool_var(f"b_{i}") for i in range(n)]
    model.add(b[0]==0)
    for i in range(n):
        ns=list(neighbors(d,i))
        model.add(sum(s*(1-2*b[j]) for j,s in ns)==root*(1-2*b[i]))
    solver=cp_model.CpSolver()
    solver.parameters.max_time_in_seconds=seconds
    solver.parameters.num_search_workers=8
    status=solver.solve(model)
    result=dict(d=d,status=solver.status_name(status),wall_time=solver.wall_time)
    if status in [cp_model.OPTIMAL,cp_model.FEASIBLE]:
        f=np.array([1-2*solver.value(x) for x in b],dtype=np.int64)
        assert np.array_equal(apply_h(f,d),root*f)
        result["witness"]=f.tolist()
    print(json.dumps(result),flush=True)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--d",type=int,default=9)
    p.add_argument("--seconds",type=float,default=60)
    args=p.parse_args()
    solve(args.d,args.seconds)
