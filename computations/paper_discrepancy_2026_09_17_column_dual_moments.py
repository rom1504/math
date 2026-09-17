"""Exact energy/covariance diagnostics of frozen actual-code minimax duals."""
from fractions import Fraction as F
from pathlib import Path
import importlib.util
import itertools
import json
import math
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("stored",ROOT/"computations/transfer_adversary_minimizer_isotropy_2026_09_06.py")
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
matrices={name:np.asarray(a,dtype=np.int64) for name,a,_ in module.cases()}
data=json.loads((ROOT/"tmp/paper_portfolio_2026_09_17/discrepancy/full_column_game_certificates.json").read_text())
reports=[]
for row in data["cases"]:
    if row["deficit_window"] or row["isotropy_required"] or row["n"]<6:
        continue
    n=row["n"]
    a=matrices[row["case"]]
    edges=list(itertools.combinations(range(n),2))
    signs=np.asarray([a[i,j] for i,j in edges])
    words=np.asarray([(1,)+h for h in itertools.product((-1,1),repeat=n-1)],dtype=np.int64)
    features=np.asarray([words[:,i]*words[:,j] for i,j in edges]).T
    energy=features@signs
    atoms=row["dual_code_law"]
    code=np.asarray([atom["word"] for atom in atoms],dtype=np.int64)
    weights=[F(atom["weight"]) for atom in atoms]
    assert sum(weights)==1
    den=math.lcm(*(w.denominator for w in weights))
    nums=np.asarray([w.numerator*(den//w.denominator) for w in weights],dtype=object)
    code_features=np.asarray([code[:,i]*code[:,j] for i,j in edges]).T
    code_energy=code_features@signs
    sectors=np.sign(code_energy)
    second=code_features.astype(object).T@nums
    signed_second=code_features.astype(object).T@(nums*sectors)
    signed_mass=F(int(nums@sectors),den)
    unsigned_energy=F(int(second@signs),den)
    signed_energy=F(int(signed_second@signs),den)
    assert signed_energy==row["cap"]
    covariance=np.eye(n)
    signed_covariance=float(signed_mass)*np.eye(n)
    signed_covariance_integer=int(nums@sectors)*np.eye(n,dtype=object)
    for k,(i,j) in enumerate(edges):
        covariance[i,j]=covariance[j,i]=float(F(int(second[k]),den))
        signed_covariance[i,j]=signed_covariance[j,i]=float(F(int(signed_second[k]),den))
        signed_covariance_integer[i,j]=signed_covariance_integer[j,i]=int(signed_second[k])
    radial_signed=F(row["cap"],len(edges))
    residual=[F(int(v),den)-radial_signed*int(s) for v,s in zip(signed_second,signs)]
    response=np.abs(words@code.T).astype(object)@nums
    minimum=int(min(response))
    own_response=np.abs(code@code.T).astype(object)@nums
    squared=signed_covariance_integer@signed_covariance_integer
    square_scalar=(str(F(int(squared[0,0]),den*den))
                   if np.array_equal(squared,int(squared[0,0])*np.eye(n,dtype=object)) else None)
    minimizing=energy[response==minimum]
    primal=row["primal_law"]
    primal_energies={}
    for atom in primal:
        word=np.asarray(atom["word"])
        value=int(sum(a[i,j]*word[i]*word[j] for i,j in edges))
        primal_energies[value]=primal_energies.get(value,F(0))+F(atom["weight"])
    reports.append({"case":row["case"],"dual_support":len(code),"value":str(F(minimum,den)),
                    "signed_mass":str(signed_mass),"signed_energy_exact":str(signed_energy),
                    "own_code_response_range_exact":[str(F(int(min(own_response)),den)),str(F(int(max(own_response)),den))],
                    "signed_covariance_square_scalar_exact":square_scalar,
                    "signed_covariance_exact":[[str(F(int(v),den)) for v in line] for line in signed_covariance_integer],
                    "unsigned_energy_exact":str(unsigned_energy),
                    "unsigned_covariance_max_offdiag":str(max(abs(F(int(v),den)) for v in second)),
                    "unsigned_covariance_spectrum_diagnostic":np.linalg.eigvalsh(covariance).tolist(),
                    "signed_covariance_spectrum_diagnostic":np.linalg.eigvalsh(signed_covariance).tolist(),
                    "signed_radial_residual_max_exact":str(max(map(abs,residual))),
                    "minimizing_column_count":len(minimizing),
                    "minimizing_column_energy_histogram":{str(int(e)):int(np.sum(minimizing==e)) for e in np.unique(minimizing)},
                    "chosen_primal_column_energy_law":{str(e):str(w) for e,w in sorted(primal_energies.items())}})
print(json.dumps({"status":"PASS exact moment identities; spectra diagnostic only","cases":reports},indent=2))
