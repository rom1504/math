"""Actual n=7 minimizers: spectral screen, then numerical SDP per invariant bin.

Invariant bins are only a candidate selection mechanism, NOT a claimed complete
switching-isomorphism classification or a certificate for omitted candidates.
"""
import json
from pathlib import Path
import numpy as np
import cvxpy as cp

n=7
pairs=[(i,j) for i in range(1,n) for j in range(i+1,n)]
spins=1-2*((np.arange(1<<(n-1))[:,None]>>np.arange(n-1))&1)
spins=np.column_stack([np.ones(len(spins),dtype=int),spins])
tests=np.array([spins[:,i]*spins[:,j] for i,j in pairs],dtype=np.int16)
star=spins[:,1:].sum(axis=1)
bins={}
caps=[]
for lo in range(0,1<<len(pairs),4096):
    ids=np.arange(lo,min(lo+4096,1<<len(pairs)))
    signs=(1-2*((ids[:,None]>>np.arange(len(pairs)))&1)).astype(np.int16)
    energies=signs@tests+star
    for idx,s,e in zip(ids,signs,energies):
        if np.max(np.abs(e))!=9:continue
        A=np.ones((n,n),dtype=np.int64)-np.eye(n,dtype=np.int64)
        for (i,j),v in zip(pairs,s):A[i,j]=A[j,i]=v
        powers=np.eye(n,dtype=np.int64)
        key=[]
        for k in range(1,n+1):
            powers=powers@A
            key.append(int(np.trace(powers)))
        key=tuple(key)
        if key not in bins:bins[key]=[int(idx),A,0]
        bins[key][2]+=1

out=[]
for key,(idx,A,count) in sorted(bins.items()):
    R=cp.Variable((n,n),symmetric=True)
    vals=[]
    for polarity in [1,-1]:
        pr=cp.Problem(cp.Maximize(polarity*cp.sum(cp.multiply(A,R))/2),[R>>0,cp.diag(R)==1])
        val=pr.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_feas=1e-9,tol_gap_rel=1e-9)
        vals.append(dict(polarity=polarity,value=float(val),R=R.value.tolist(),status=pr.status))
    row=dict(id=idx,power_traces=key,count=count,SDP=max(v['value'] for v in vals),results=vals)
    out.append(row)
    print(json.dumps({k:v for k,v in row.items() if k!='results'}),flush=True)
Path(__file__).with_suffix('.json').write_text(json.dumps(dict(total=sum(v[2] for v in bins.values()),bins=out),indent=2)+'\n')
