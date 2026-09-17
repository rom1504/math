"""Recover the remaining n12 radial-mixture endpoint from its active face."""
from fractions import Fraction as F
from pathlib import Path
import importlib.util
import itertools
import json
import math
import numpy as np
import sympy as sp

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("stored",ROOT/"computations/transfer_adversary_minimizer_isotropy_2026_09_06.py")
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
matrix=next(a for _,a,_ in module.cases() if len(a)==12)
data=json.loads((ROOT/"tmp/paper_portfolio_2026_09_17/discrepancy/radial_correlation_mixture_certificates.json").read_text())
row=next(a for a in data["cases"] if a["n"]==12)
mixture=row["isotropic_mixture"]
n=12
words=np.asarray([(1,)+h for h in itertools.product((-1,1),repeat=n-1)],dtype=np.int64)
edges=list(itertools.combinations(range(n),2))
features=np.asarray([words[:,i]*words[:,j] for i,j in edges]).T
block=np.column_stack((np.ones(len(words),dtype=np.int64),features))
code=np.asarray([a["word"] for a in mixture["optimized_dual_ground_law"]],dtype=np.int64)
overlap=np.abs(words@code.T)
weight=F(mixture["positive_radial_law_weight_exact"])
assert weight==F(1,2)
count=block.shape[1]
# Multiply both query inequalities by two to keep all coefficients integral.
coefficients=np.zeros((2*len(words),2*count+len(code)),dtype=np.int64)
coefficients[:len(words),:count]=-2*block
coefficients[len(words):,count:2*count]=-2*block
coefficients[:len(words),2*count:]=overlap
coefficients[len(words):,2*count:]=overlap
approximate=np.asarray([float(F(v)) for v in mixture["optimized_dual_equality_coefficients"]]+
                       [float(F(a["weight"])) for a in mixture["optimized_dual_ground_law"]])
tight=np.flatnonzero(np.abs(coefficients@approximate)<1e-7)
rows=[list(map(int,coefficients[i]))+[0] for i in tight]
rows.append([0]*(2*count)+[1]*len(code)+[1])
signs=[matrix[i][j] for i,j in edges]
rhs=[F(1)]+[F(a,4) for a in signs]+[F(1)]+[F(-a,4) for a in signs]
upper=F(mixture["optimized_ground_response_upper_exact"])
rows.append([sp.Rational(v.numerator,v.denominator) for v in rhs]+[0]*len(code)+[sp.Rational(upper.numerator,upper.denominator)])
print("face",len(rows),len(rows[0]),flush=True)
reduced_domain,pivots=sp.polys.matrices.DomainMatrix.from_Matrix(sp.Matrix(rows)).rref()
reduced=reduced_domain.to_Matrix()
columns=len(approximate)
assert columns not in pivots
free=[j for j in range(columns) if j not in pivots]
values=[F(0)]*columns
for j in free:
    values[j]=F(float(approximate[j])).limit_denominator(1000000)
for i,j in enumerate(pivots):
    values[j]=F(reduced[i,columns])-sum(F(reduced[i,k])*values[k] for k in free)
den=math.lcm(*(v.denominator for v in values))
nums=np.asarray([v.numerator*(den//v.denominator) for v in values],dtype=object)
margins=coefficients.astype(object)@nums
pi=values[2*count:]
objective=sum(a*b for a,b in zip(rhs,values))
print("result","free",len(free),"minpi",str(min(pi)),"sum",str(sum(pi)),
      "margin",str(F(int(min(margins)),den)),"objective",str(objective),flush=True)
assert min(pi)>=0 and sum(pi)<=1 and min(margins)>=0 and objective==upper
print(json.dumps({"case":row["case"],"optimized_ground_response_exact":str(upper),
                  "optimized_dual_equality_coefficients":[str(v) for v in values[:2*count]],
                  "optimized_dual_ground_law":[{"word":word.tolist(),"weight":str(w)} for word,w in zip(code,pi) if w]},separators=(",",":")))
