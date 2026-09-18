"""Numerical MILP diagnostic: exact optimal small children, bridge free.

Solver status is numerical evidence, not a formal lower certificate.
All attained caps are replayed by exhaustive integer spin evaluation.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import Bounds, LinearConstraint, milp
from scipy.sparse import csc_matrix

out = []
seeds = [(4,1),(5,13),(6,220),(7,691)]
for n, idx in seeds:
    child=np.ones((n,n),dtype=np.int16)-np.eye(n,dtype=np.int16)
    for k,(i,j) in enumerate((i,j) for i in range(1,n) for j in range(i+1,n)):
        child[i,j]=child[j,i]=1-2*((idx>>k)&1)
    N=2*n
    spins=1-2*((np.arange(1<<(N-1))[:,None]>>np.arange(N-1))&1)
    spins=np.column_stack([np.ones(len(spins),dtype=int),spins])
    bridge=np.array([spins[:,i]*spins[:,n+j] for i in range(n) for j in range(n)],dtype=float).T
    for polarity in [1,-1]:
        base=(np.einsum('bi,ij,bj->b',spins[:,:n],child,spins[:,:n])+
              polarity*np.einsum('bi,ij,bj->b',spins[:,n:],child,spins[:,n:]))//2
        offset=base-bridge.sum(axis=1)
        rows=np.vstack([np.column_stack([2*bridge,-np.ones(len(spins))]),
                        np.column_stack([-2*bridge,-np.ones(len(spins))])])
        rhs=np.r_[-offset,offset]
        c=np.r_[np.zeros(n*n),1.]
        result=milp(c,integrality=np.r_[np.ones(n*n),0],bounds=Bounds(np.zeros(n*n+1),np.r_[np.ones(n*n),np.inf]),
                    constraints=LinearConstraint(csc_matrix(rows),-np.inf,rhs),
                    options={'time_limit':45,'mip_rel_gap':0})
        row=dict(child_n=n,seed=idx,polarity=polarity,status=int(result.status),message=result.message)
        if result.x is not None:
            signs=2*np.rint(result.x[:n*n]).astype(int)-1
            cap=int(np.max(np.abs(base+bridge@signs)))
            row.update(cap=cap,dual_bound=float(result.mip_dual_bound),gap=float(result.mip_gap),bridge=signs.reshape(n,n).tolist())
        out.append(row)
        print(json.dumps(row),flush=True)
        Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
