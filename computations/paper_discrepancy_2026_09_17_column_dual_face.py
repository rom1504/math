"""Exact affine-face recovery of the two remaining isotropic shell duals."""
from fractions import Fraction as F
from pathlib import Path
import itertools
import json
import math
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
data=json.loads((ROOT/"tmp/paper_portfolio_2026_09_17/discrepancy/full_column_game_certificates.json").read_text())
for row in data["cases"]:
    if not (row["n"]==12 and row["isotropy_required"] and row["deficit_window"] in (2,4)):
        continue
    n=row["n"]
    words=np.asarray([(1,)+h for h in itertools.product((-1,1),repeat=n-1)],dtype=np.int64)
    edges=list(itertools.combinations(range(n),2))
    features=np.asarray([words[:,i]*words[:,j] for i,j in edges]).T
    code=np.asarray([a["word"] for a in row["dual_code_law"]],dtype=np.int64)
    overlap=np.abs(words@code.T)
    approximate=np.asarray([float(F(v)) for v in row["dual_equality_coefficients"][1:]]+
                           [float(F(a["weight"])) for a in row["dual_code_law"]])
    coefficients=np.column_stack((-features,overlap))
    target=F(row["upper_exact"])
    residual=coefficients@approximate-float(target)
    tight=np.flatnonzero(np.abs(residual)<1e-7)
    matrix=[list(map(int,coefficients[i]))+[sp.Rational(target.numerator,target.denominator)] for i in tight]
    matrix.append([0]*len(edges)+[1]*len(code)+[1])
    print("recover",row["deficit_window"],"shape",len(matrix),len(matrix[0]),flush=True)
    reduced_domain,pivots=sp.polys.matrices.DomainMatrix.from_Matrix(sp.Matrix(matrix)).rref()
    reduced=reduced_domain.to_Matrix()
    columns=len(approximate)
    assert columns not in pivots
    free=[j for j in range(columns) if j not in pivots]
    values=[F(0)]*columns
    for j in free:
        values[j]=F(float(approximate[j])).limit_denominator(1000000)
    for i,j in enumerate(pivots):
        values[j]=F(reduced[i,columns])-sum(F(reduced[i,k])*values[k] for k in free)
    den=math.lcm(*(v.denominator for v in values+[target]))
    nums=np.asarray([v.numerator*(den//v.denominator) for v in values],dtype=object)
    margins=coefficients.astype(object)@nums-target.numerator*(den//target.denominator)
    pi=values[len(edges):]
    print("result",row["deficit_window"],"free",len(free),"minpi",str(min(pi)),
          "sum",str(sum(pi)),"margin",str(F(int(min(margins)),den)),flush=True)
    if min(pi)>=0 and sum(pi)<=1 and min(margins)>=0:
        print(json.dumps({"case":row["case"],"deficit_window":row["deficit_window"],
                          "lower_exact":str(target),"dual_equality_coefficients":[str(target)]+[str(v) for v in values[:len(edges)]],
                          "dual_code_law":[{"word":word.tolist(),"weight":str(w)} for word,w in zip(code,pi) if w]},separators=(",",":")),flush=True)
