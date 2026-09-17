"""Finite mode-law LP on exact stored ground codes.

The LP is numerical discovery; rationalized primal/dual inequalities are
checked exactly where reconstruction succeeds. No asymptotic claim is
inferred from finite seed data. Uses tracked source witnesses only.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import importlib.util
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linprog
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]


def rational_vector(values, max_denominator=1000000):
    return [Fraction(float(v)).limit_denominator(max_denominator) for v in values]


def solve_case(name, matrix, provenance, slack):
    a = np.asarray(matrix, dtype=np.int64)
    k = len(a)
    z = np.asarray([(1,) + v for v in itertools.product((-1, 1), repeat=k-1)], dtype=np.int64)
    energy = np.einsum("bi,ij,bj->b", z, a, z)//2
    cap = int(np.abs(energy).max())
    assert cap == provenance["minimum_cap_imported"]
    code = z[np.abs(energy) >= cap-slack]
    pair = np.asarray([z[:, i]*z[:, j] for i,j in itertools.combinations(range(k),2)])
    equality = np.vstack((np.ones(len(z), dtype=np.int64), pair))
    costs = np.abs(code @ z.T)
    ae = np.column_stack((equality, np.zeros(len(equality))))
    be = np.zeros(len(equality)); be[0] = 1
    au = np.column_stack((costs, -np.ones(len(code))))
    objective = np.zeros(len(z)+1); objective[-1] = 1
    solution = linprog(objective, A_eq=ae, b_eq=be, A_ub=au,
                       b_ub=np.zeros(len(code)), bounds=(0,None), method="highs")
    if not solution.success:
        raise RuntimeError(solution.message)
    support = np.flatnonzero(solution.x[:-1] > 1e-9)
    weights = rational_vector(solution.x[support])
    mean = [sum(w*int(v) for w,v in zip(weights,row)) for row in equality[:,support]]
    primal_ok = mean == [Fraction(1)] + [Fraction(0)]*(len(equality)-1)
    repair_status = "not needed"
    if not primal_ok:
        # The numerical vertex's active equations have integer coefficients.
        # Solve those equations over Q instead of claiming that rounded
        # floating weights are exactly isotropic.
        active=np.flatnonzero(np.abs(costs @ solution.x[:-1]-solution.fun)<1e-7)
        rows=np.vstack((np.column_stack((equality[:,support],np.zeros(len(equality),dtype=int))),
                        np.column_stack((costs[active][:,support],-np.ones(len(active),dtype=int)))))
        target=[1]+[0]*(len(rows)-1)
        try:
            answer,parameters=sp.Matrix(rows.tolist()).gauss_jordan_solve(sp.Matrix(target))
            if len(parameters):
                repair_status="active system has free parameters; not certified"
            else:
                weights=[Fraction(int(v.p),int(v.q)) for v in answer[:-1]]
                mean=[sum(w*int(v) for w,v in zip(weights,row)) for row in equality[:,support]]
                primal_ok=(min(weights)>=0 and mean==[Fraction(1)]+[Fraction(0)]*(len(equality)-1))
                repair_status="exact active-system reconstruction" if primal_ok else "reconstruction infeasible"
        except (ValueError, ZeroDivisionError):
            repair_status="active-system reconstruction failed; no exact claim"
    primal_bound = max(sum(w*int(v) for w,v in zip(weights,row)) for row in costs[:,support])
    dual_eq = rational_vector(solution.eqlin.marginals)
    dual_ub = [min(Fraction(0),v) for v in rational_vector(solution.ineqlin.marginals)]
    dual_scale=max(Fraction(1),-sum(dual_ub))
    dual_eq=[v/dual_scale for v in dual_eq]
    dual_ub=[v/dual_scale for v in dual_ub]
    # Dual: eq^T a + cost^T u <=0; -sum u<=1; u<=0.
    dual_ok = all(v<=0 for v in dual_ub) and -sum(dual_ub)<=1
    maximum_violation=Fraction(0)
    if dual_ok:
        for index in range(len(z)):
            value = sum(w*int(v) for w,v in zip(dual_eq,equality[:,index]))
            value += sum(w*int(v) for w,v in zip(dual_ub,costs[:,index]))
            maximum_violation=max(maximum_violation,value)
        # The first equality row is identically one, so decreasing its
        # dual coefficient repairs every residual exactly.
        dual_eq[0]-=maximum_violation
    return {"case":name,"k":k,"cap":cap,"slack":slack,
            "code_size_projective":len(code),"lp_value":float(solution.fun),
            "normalized_mean_coefficient":float(solution.fun)/np.sqrt(k),
            "uniform_sign_mean":float(costs[0].mean()),
            "primal_rational_verified":primal_ok,
            "primal_repair":repair_status,
            "primal_upper_bound":str(primal_bound) if primal_ok else None,
            "dual_rational_verified":dual_ok,
            "dual_repair_decrease":str(maximum_violation),
            "dual_lower_bound":str(dual_eq[0]) if dual_ok else None,
            "support":[{"mode":z[i].tolist(),"weight":str(w)} for i,w in zip(support,weights)],
            "dual_equalities":[str(w) for w in dual_eq],
            "dual_inequalities":[str(w) for w in dual_ub],
            "scope":"finite moment-law certificate; original optimality imported, cap replayed"}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-order",type=int,default=14)
    parser.add_argument("--slack",type=int,default=0)
    parser.add_argument("--output",type=Path)
    args=parser.parse_args()
    source=ROOT/"computations/transfer_adversary_minimizer_isotropy_2026_09_06.py"
    spec=importlib.util.spec_from_file_location("stored_mode_inputs",source)
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    reports=[]
    for name,matrix,provenance in module.cases():
        if len(matrix)%2 or len(matrix)>args.max_order: continue
        report=solve_case(name,matrix,provenance,args.slack)
        reports.append(report)
        print(json.dumps({key:value for key,value in report.items()
                          if key not in ("support","dual_equalities","dual_inequalities")}),flush=True)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(json.dumps({"cases":reports},indent=2)+"\n")


if __name__=="__main__":
    main()
