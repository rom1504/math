#!/usr/bin/env python3
"""Bounded CP-SAT feasibility for all chiral full signings of order 16.

The copy relation B~A is deliberately omitted. Pair rotations normalize
d=+1; paired switches normalize A[0,i]=+1. Permutations sort the first
row of B, and exchanging halves complements its signs after re-gauging.
Any INFEASIBLE result is solver-certified only, without a proof object.
"""

import argparse
import json
from pathlib import Path
import time

import numpy as np
from ortools.sat.python import cp_model


ROOT = Path(__file__).resolve().parents[1]


def projective_spins(n):
    ids = np.arange(1 << (n-1), dtype=np.int64)
    return np.column_stack((np.ones(len(ids), dtype=np.int64),
                            1 - 2 * ((ids[:, None] >> np.arange(n-1)) & 1)))


def normalize(a, b, d):
    n = len(a)
    c = b + np.diag(d)
    matrix = np.block([[a, c], [c, -a]])
    transform = np.eye(2*n, dtype=np.int64)
    for i in range(n):
        if d[i] < 0:
            transform[np.ix_([i, i+n], [i, i+n])] = np.array([[0, -1], [1, 0]])
    matrix = transform.T @ matrix @ transform
    aa = matrix[:n, :n]
    bb = matrix[:n, n:] - np.eye(n, dtype=np.int64)
    signs = np.r_[1, aa[0, 1:]]
    aa = aa * signs[:, None] * signs[None, :]
    bb = bb * signs[:, None] * signs[None, :]
    if np.count_nonzero(bb[0, 1:] < 0) > (n-1)//2:
        aa = -aa
        signs = np.r_[1, -np.ones(n-1, dtype=np.int64)]
        aa = aa * signs[:, None] * signs[None, :]
        bb = bb * signs[:, None] * signs[None, :]
    p = np.r_[0, 1 + np.argsort(bb[0, 1:], kind="stable")]
    aa, bb = aa[np.ix_(p, p)], bb[np.ix_(p, p)]
    assert np.all(aa[0, 1:] == 1)
    assert np.all(np.diff(bb[0, 1:]) >= 0)
    assert np.count_nonzero(bb[0, 1:] < 0) <= (n-1)//2
    return aa, bb


def build_constraints(n, cap):
    xs = projective_spins(n)
    variables = [("A", i, j) for i in range(1, n) for j in range(i+1, n)]
    variables += [("B", i, j) for i in range(n) for j in range(i+1, n)]
    constraints = {}
    for x in xs:
        for y in xs:
            coeff = np.array([(x[i]*x[j]-y[i]*y[j])//2 if kind == "A"
                              else (x[i]*y[j]+x[j]*y[i])//2
                              for kind, i, j in variables], dtype=np.int64)
            fixed = sum((x[0]*x[j]-y[0]*y[j])//2 for j in range(1, n))
            const = fixed + int(x @ y)//2 - int(coeff.sum())
            # H_D/2 = const + 2 * sum(coeff_i * bit_i), bit_i=(edge_i+1)/2.
            lower = (-cap//2 - const + 1)//2
            upper = (cap//2 - const)//2
            nz = np.flatnonzero(coeff)
            if len(nz) and coeff[nz[0]] < 0:
                coeff = -coeff
                lower, upper = -upper, -lower
            key = tuple(map(int, coeff))
            if key in constraints:
                old_l, old_u = constraints[key]
                constraints[key] = (max(lower, old_l), min(upper, old_u))
            else:
                constraints[key] = (lower, upper)
    rows = np.array(list(constraints), dtype=np.int8)
    bounds = np.array(list(constraints.values()), dtype=np.int64)
    return variables, rows, bounds


def exact_cap(a, b):
    x = projective_spins(len(a))
    q = np.einsum("bi,ij,bj->b", x, a, x)//2
    values = q[:, None] - q[None, :] + x @ (b + np.eye(len(a), dtype=np.int64)) @ x.T
    return int(abs(values).max())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cap", type=int, default=28)
    parser.add_argument("--seconds", type=float, default=600)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--seed", type=int, default=20260918)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--hint", type=Path, default=ROOT / "computations/results/twisted_chiral_2026_09_18_profile_width_m8_class0_checkpoint.json")
    args = parser.parse_args()
    started = time.time()
    n = 8
    variables, rows, bounds = build_constraints(n, args.cap)
    model = cp_model.CpModel()
    bits = [model.NewBoolVar("%s_%d_%d" % item) for item in variables]
    for row, (lower, upper) in zip(rows, bounds):
        model.AddLinearConstraint(sum(int(c)*bits[i] for i, c in enumerate(row) if c), int(lower), int(upper))
    broot = [bits[variables.index(("B", 0, j))] for j in range(1, n)]
    for left, right in zip(broot, broot[1:]):
        model.Add(left <= right)
    model.Add(sum(broot) >= (n-1+1)//2)
    hint_cap = None
    if args.hint.exists():
        hint = json.loads(args.hint.read_text())
        ha, hb = normalize(np.asarray(hint["child_matrix"]),
                           np.asarray(hint["bridge_hollow_matrix"]), np.asarray(hint["d"]))
        hint_cap = exact_cap(ha, hb)
        assert hint_cap == hint["cap"]
        for variable, (kind, i, j) in zip(bits, variables):
            edge = (ha if kind == "A" else hb)[i, j]
            model.AddHint(variable, int((edge+1)//2))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(args.output.with_suffix(".npz"), rows=rows, bounds=bounds,
                        variable_kind=np.array([v[0] for v in variables]),
                        variable_i=np.array([v[1] for v in variables]),
                        variable_j=np.array([v[2] for v in variables]))
    metadata = {"order": 16, "cap_target": args.cap, "variables": len(bits),
                "distinct_cap_rows": len(rows), "absolute_energy_orbit_reps": 16384,
                "gauge": "d=+1, A0i=+1, sorted B0i, at most three negative B0i",
                "scope": "all symmetric chiral doubles with arbitrary independent A,B; no B~A restriction",
                "seconds": args.seconds, "workers": args.workers, "seed": args.seed,
                "hint_cap": hint_cap, "solver_proof_object": False,
                "status": "RUNNING", "model_stats": model.ModelStats()}
    args.output.write_text(json.dumps(metadata, indent=2)+"\n")
    print(json.dumps(metadata), flush=True)
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = args.seconds
    solver.parameters.num_search_workers = args.workers
    solver.parameters.random_seed = args.seed
    solver.parameters.log_search_progress = True
    status = solver.Solve(model)
    metadata.update(status=solver.StatusName(status), elapsed_seconds=time.time()-started,
                    response_stats=solver.ResponseStats())
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        a = np.ones((n, n), dtype=np.int64)-np.eye(n, dtype=np.int64)
        b = np.zeros((n, n), dtype=np.int64)
        for variable, (kind, i, j) in zip(bits, variables):
            matrix = a if kind == "A" else b
            matrix[i, j] = matrix[j, i] = 2*solver.Value(variable)-1
        cap = exact_cap(a, b)
        assert cap <= args.cap
        metadata.update(A=a.tolist(), B=b.tolist(), d=[1]*n, exact_witness_cap=cap)
    args.output.write_text(json.dumps(metadata, indent=2)+"\n")
    print(json.dumps(metadata), flush=True)


if __name__ == "__main__":
    main()
