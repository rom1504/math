"""Full physical-column response games on exact finite near-level codes.

Both unrestricted and isotropic laws are included. Every primal and dual
is replayed rationally, with certified rational intervals if exact
endpoint reconstruction fails. The matrices' optimum labels are imported.
"""

from __future__ import annotations

import argparse
from fractions import Fraction as F
import importlib.util
import itertools
import json
import math
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.optimize import linprog


ROOT = Path(__file__).resolve().parents[1]


def load_module(name,path):
    spec = importlib.util.spec_from_file_location(name,ROOT/path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def integerize(values):
    den = math.lcm(*(value.denominator for value in values))
    return den,np.asarray([value.numerator*(den//value.denominator) for value in values],dtype=object)


def recover_isotropic_dual_face(result,eq,overlap,upper):
    """Recover a rational endpoint on the complete numerical active face."""
    support=np.flatnonzero(-result.ineqlin.marginals>1e-9)
    coefficients=np.column_stack((-eq[1:].T,overlap[support].T))
    approximate=np.r_[result.eqlin.marginals[1:],-result.ineqlin.marginals[support]]
    tight=np.flatnonzero(np.abs(coefficients@approximate-float(upper))<1e-7)
    target=sp.Rational(upper.numerator,upper.denominator)
    rows=[list(map(int,coefficients[i]))+[target] for i in tight]
    rows.append([0]*(len(eq)-1)+[1]*len(support)+[1])
    reduced_domain,pivots=sp.polys.matrices.DomainMatrix.from_Matrix(sp.Matrix(rows)).rref()
    columns=len(approximate)
    if columns in pivots:
        return None
    reduced=reduced_domain.to_Matrix()
    free=[j for j in range(columns) if j not in pivots]
    values=[F(0)]*columns
    for j in free:
        values[j]=F(float(approximate[j])).limit_denominator(1000000)
    for i,j in enumerate(pivots):
        values[j]=F(reduced[i,columns])-sum(F(reduced[i,k])*values[k] for k in free)
    pi=[F(0)]*len(overlap)
    for i,v in zip(support,values[len(eq)-1:]):
        pi[i]=v
    dual=[upper]+values[:len(eq)-1]
    den,nums=integerize(dual+pi)
    margins=overlap.T.astype(object)@nums[len(eq):]-eq.T.astype(object)@nums[:len(eq)]
    if min(pi)>=0 and sum(pi)<=1 and min(margins)>=0:
        return dual,pi
    return None


def certify(result,eq,overlap,words,code,module,full,unrestricted=None):
    indices = np.flatnonzero(result.x[:-1] > 1e-8)
    weights = None
    for bound in (1000,1000000,1000000000):
        proposed = [F(float(result.x[i])).limit_denominator(bound) for i in indices]
        den,nums = integerize(proposed)
        actual = eq[:,indices].astype(object)@nums
        if int(actual[0]) == den and all(int(v)==0 for v in actual[1:]):
            weights = proposed
            break
    if weights is None:
        rows = [row.tolist()+[0] for row in eq[:,indices]]
        rhs = [1]+[0]*(len(eq)-1)
        tight = np.flatnonzero(np.abs(overlap@result.x[:-1]-result.x[-1]) < 1e-7)
        rows += [overlap[i,indices].tolist()+[-1] for i in tight]
        rhs += [0]*len(tight)
        values = module.exact_linear_solution(rows,rhs,np.r_[result.x[indices],result.x[-1]])
        assert values is not None and min(values)>=0
        weights = values[:-1]
    den,nums = integerize(weights)
    response = overlap[:,indices].astype(object)@nums
    upper = F(int(max(response)),den)
    exact = None
    for bound in (1000,1000000,1000000000):
        dual = [F(float(v)).limit_denominator(bound) if abs(v)>1e-9 else F(0)
                for v in result.eqlin.marginals]
        pi = [F(float(-v)).limit_denominator(bound) if abs(v)>1e-9 else F(0)
              for v in result.ineqlin.marginals]
        if min(pi)<0 or sum(pi)>1 or dual[0]!=upper:
            continue
        den,nums = integerize(dual+pi)
        margins = overlap.T.astype(object)@nums[len(eq):]-eq.T.astype(object)@nums[:len(eq)]
        if min(margins)>=0:
            exact = (dual,pi)
            break
    status = "exact rational primal and all-column dual agree"
    if (exact is None and len(eq)>1 and unrestricted is not None
        and F(unrestricted["lower_exact"])==upper):
        code_weights={tuple(a["word"]):F(a["weight"]) for a in unrestricted["dual_code_law"]}
        exact=([upper]+[F(0)]*(len(eq)-1),[code_weights.get(tuple(word),F(0)) for word in code])
        status="exact rational primal and transferred unrestricted all-column dual agree"
    if exact is None and len(eq)>1:
        exact=recover_isotropic_dual_face(result,eq,overlap,upper)
        if exact is not None:
            status="exact rational primal and exact affine-face all-column dual agree; direct rounding failed"
    if exact is None and len(eq)==1:
        code_index = {tuple(word):i for i,word in enumerate(code)}
        if all(tuple(words[i]) in code_index for i in indices):
            pi = [F(0)]*len(code)
            for i,w in zip(indices,weights):
                pi[code_index[tuple(words[i])]] += w
            den,nums = integerize([upper]+pi)
            margins = overlap.T.astype(object)@nums[1:]-nums[0]
            if min(margins)>=0:
                exact = ([upper],pi)
                status = "exact rational primal is also an all-column code-side dual"
    if exact is None:
        precision = 10**12
        dual = [F(round(float(v)*precision),precision) for v in result.eqlin.marginals]
        pi = [F(max(0,round(float(-v)*precision)),precision) for v in result.ineqlin.marginals]
        scale = max(F(1),sum(pi))
        dual = [v/scale for v in dual]
        pi = [v/scale for v in pi]
        den,nums = integerize(dual+pi)
        margins = overlap.T.astype(object)@nums[len(eq):]-eq.T.astype(object)@nums[:len(eq)]
        dual[0] += min(F(0),F(int(min(margins)),den))
        exact = (dual,pi)
        status = "exact rational interval after failed endpoint recovery and pointwise dual correction"
    dual,pi = exact
    den,nums = integerize(dual+pi)
    margins = overlap.T.astype(object)@nums[len(eq):]-eq.T.astype(object)@nums[:len(eq)]
    assert min(margins)>=0 and min(pi)>=0 and sum(pi)<=1
    lower = dual[0]
    assert lower <= upper and upper-lower < F(1,1000000)
    report = {"lower_exact":str(lower),"upper_exact":str(upper),
              "normalized_upper":float(upper)/math.sqrt(words.shape[1]),
              "status":status,"primal_support_size":len(indices),
              "dual_code_support_size":sum(v>0 for v in pi)}
    if full:
        report["primal_law"] = [{"word":words[i].tolist(),"weight":str(w)} for i,w in zip(indices,weights)]
        report["dual_equality_coefficients"] = [str(v) for v in dual]
        report["dual_code_law"] = [{"word":word.tolist(),"weight":str(w)} for word,w in zip(code,pi) if w]
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--orders",default="8,10")
    parser.add_argument("--windows",default="0,2,4")
    parser.add_argument("--full",action="store_true")
    parser.add_argument("--compact",action="store_true")
    args = parser.parse_args()
    orders = {int(value) for value in args.orders.split(",")}
    windows = [int(value) for value in args.windows.split(",")]
    module = load_module("stored_minimizers","computations/transfer_adversary_minimizer_isotropy_2026_09_06.py")
    reports = []
    for name,matrix,provenance in module.cases():
        n = len(matrix)
        if n not in orders:
            continue
        a = np.asarray(matrix,dtype=np.int64)
        words = np.asarray([(1,)+rest for rest in itertools.product((-1,1),repeat=n-1)],dtype=np.int64)
        edges = list(itertools.combinations(range(n),2))
        features = np.asarray([words[:,i]*words[:,j] for i,j in edges],dtype=np.int64).T
        energies = features@np.asarray([a[i,j] for i,j in edges])
        cap = int(np.abs(energies).max())
        assert cap == provenance["minimum_cap_imported"]
        for window in windows:
            code = words[cap-np.abs(energies) <= window]
            overlap = np.abs(code@words.T)
            unrestricted = None
            for isotropic in (False,True):
                eq = (np.vstack((np.ones(len(words),dtype=np.int64),features.T))
                      if isotropic else np.ones((1,len(words)),dtype=np.int64))
                result = linprog(np.r_[np.zeros(len(words)),1.0],
                                 A_eq=np.column_stack((eq,np.zeros(len(eq)))),
                                 b_eq=np.r_[1.0,np.zeros(len(eq)-1)],
                                 A_ub=np.column_stack((overlap,-np.ones(len(code)))),
                                 b_ub=np.zeros(len(code)),bounds=(0,None),method="highs")
                assert result.success,result.message
                report = {"case":name,"n":n,"cap":cap,"deficit_window":window,
                          "projective_code_size":len(code),"isotropy_required":isotropic,
                          "numerical_optimum":float(result.x[-1])}
                report.update(certify(result,eq,overlap,words,code,module,True,unrestricted))
                if (isotropic and unrestricted is not None
                    and F(unrestricted["lower_exact"]) == F(report["upper_exact"])):
                    report["lower_exact"] = unrestricted["lower_exact"]
                    report["dual_equality_coefficients"] = (
                        unrestricted["dual_equality_coefficients"]+["0"]*len(edges))
                    report["dual_code_law"] = unrestricted["dual_code_law"]
                    report["dual_code_support_size"] = unrestricted["dual_code_support_size"]
                    report["endpoint_recovery_before_transfer"] = report["status"]
                    report["status"] = "exact rational primal and transferred unrestricted all-column dual agree"
                if not isotropic:
                    unrestricted = dict(report)
                if not args.full:
                    for key in ("primal_law","dual_equality_coefficients","dual_code_law"):
                        report.pop(key)
                reports.append(report)
    print(json.dumps({"status":"PASS full physical-column games",
                      "global_minimizer_provenance_imported":True,"cases":reports},
                     indent=None if args.compact else 2,
                     separators=(",", ":") if args.compact else None))


if __name__ == "__main__":
    main()
