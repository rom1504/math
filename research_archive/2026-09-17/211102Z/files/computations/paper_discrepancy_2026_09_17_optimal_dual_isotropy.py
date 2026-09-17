"""Can an unrestricted response-game maximizing query law be isotropic?"""
from fractions import Fraction as F
from pathlib import Path
import importlib.util
import itertools
import json
import math
import numpy as np
from scipy.optimize import linprog

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("stored",ROOT/"computations/transfer_adversary_minimizer_isotropy_2026_09_06.py")
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
data=json.loads((ROOT/"tmp/paper_portfolio_2026_09_17/discrepancy/full_column_game_certificates.json").read_text())
reports=[]
for name,matrix,provenance in module.cases():
    n=len(matrix)
    row=next((r for r in data["cases"] if r["case"]==name and r["deficit_window"]==0 and not r["isotropy_required"]),None)
    if row is None:
        continue
    a=np.asarray(matrix,dtype=np.int64)
    words=np.asarray([(1,)+s for s in itertools.product((-1,1),repeat=n-1)],dtype=np.int64)
    edges=list(itertools.combinations(range(n),2))
    features=np.asarray([words[:,i]*words[:,j] for i,j in edges]).T
    energies=features@np.asarray([a[i,j] for i,j in edges])
    selected=np.abs(energies)==row["cap"]
    code=words[selected]
    eq=np.vstack((np.ones(len(code),dtype=np.int64),features[selected].T))
    overlap=np.abs(words@code.T)
    value=F(row["upper_exact"])
    result=linprog(np.zeros(len(code)),A_eq=eq,b_eq=np.r_[1,np.zeros(len(edges))],
                   A_ub=-overlap,b_ub=-float(value)*np.ones(len(words)),bounds=(0,None),method="highs")
    report={"case":name,"isotropic_dual_feasible_float":result.success}
    if result.success:
        support=np.flatnonzero(result.x>1e-9)
        weights=[F(float(result.x[i])).limit_denominator(1000000) for i in support]
        den=math.lcm(*(w.denominator for w in weights))
        nums=np.asarray([w.numerator*(den//w.denominator) for w in weights],dtype=object)
        second=eq[:,support].astype(object)@nums
        responses=overlap[:,support].astype(object)@nums
        assert second[0]==den and all(v==0 for v in second[1:])
        assert min(nums)>=0 and F(int(min(responses)),den)>=value
        report.update({"status":"exact isotropic optimal-query-law certificate",
                       "response_lower_exact":str(F(int(min(responses)),den)),
                       "law":[{"word":code[i].tolist(),"weight":str(w)} for i,w in zip(support,weights)]})
    else:
        report["status"]="float infeasibility only; not an exact nonexistence certificate"
    reports.append(report)
print(json.dumps({"cases":reports},indent=2))
