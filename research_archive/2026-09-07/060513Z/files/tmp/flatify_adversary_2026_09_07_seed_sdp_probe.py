"""Numerical SDP/actual-cap screen for a specific MUB-lift obstruction."""
import json
from pathlib import Path
import numpy as np
import cvxpy as cp

out=[]
for n in range(3,15):
    path=Path(f'computations/results/exact_m{n}.json')
    if not path.exists():continue
    A=np.array(json.loads(path.read_text())['matrix'],dtype=float)
    spins=1-2*((np.arange(1<<(n-1))[:,None]>>np.arange(n-1))&1)
    spins=np.column_stack([np.ones(len(spins)),spins])
    Q=float(np.max(np.abs(np.einsum('bi,ij,bj->b',spins,A,spins)))/2)
    R=cp.Variable((n,n),symmetric=True)
    results=[]
    for polarity in [1,-1]:
        prob=cp.Problem(cp.Maximize(polarity*cp.sum(cp.multiply(A,R))/2),[R>>0,cp.diag(R)==1])
        value=prob.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_feas=1e-9,tol_gap_rel=1e-9)
        results.append(dict(polarity=polarity,SDP=float(value),status=prob.status,R=R.value.tolist()))
    row=dict(n=n,source=str(path),actual_Q=Q,SDP_abs=max(x['SDP'] for x in results),
             ratio=max(x['SDP'] for x in results)/Q,threshold=float(np.pi/2),results=results)
    out.append(row)
    print(json.dumps({k:v for k,v in row.items() if k!='results'}),flush=True)
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
