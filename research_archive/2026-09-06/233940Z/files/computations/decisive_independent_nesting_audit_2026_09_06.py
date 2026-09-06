"""Replay fixed stored optimizer deletion profiles and exhaustive row extensions."""
import hashlib
import json
from pathlib import Path
import numpy as np

KNOWN={9:12,10:13,11:17,12:18,13:20,14:21}
SOURCES=[
 'computations/results/exact_m10.json',
 'computations/results/nested_10_in_11_cap17.json',
 'computations/results/heuristic_m11.json',
 'computations/results/extension_nested_m11_to_12.json',
 'computations/results/heuristic_m12.json',
 'computations/results/bridge_6_7_sign1_cap20.json',
 'computations/results/conference_completion_m13.json',
]
CACHE={}
def spins(n):
    if n not in CACHE:
        x=np.ones((1<<(n-1),n),dtype=np.int16)
        for i in range(1,n):x[:,i]=1-2*((np.arange(len(x))>>(i-1))&1)
        CACHE[n]=x
    return CACHE[n]
def energy(a):
    x=spins(len(a))
    return np.einsum('bi,ij,bj->b',x,a,x,dtype=np.int64)//2
def cap(a):return int(abs(energy(a)).max())
def extension(a):
    x=spins(len(a)); q=abs(energy(a)); best=10**9; row=None; count=0
    for start in range(0,len(x),128):
        vals=np.max(q[None,:]+abs(x[start:start+128]@x.T),axis=1)
        v=int(vals.min())
        if v<best:best=v;row=x[start+int(vals.argmin())].tolist();count=0
        if v==best:count+=int((vals==best).sum())
    return best,row,count

cases=[]
for source in SOURCES:
    p=Path(source);data=json.loads(p.read_text())
    key=next(k for k in ['matrix','parent_matrix','conference_matrix'] if k in data)
    a=np.array(data[key],dtype=np.int16)
    matrices=[(source,a)]
    if 'nested_10_in_11' in source:matrices.append((source+' [first10 child]',a[:10,:10]))
    for label,b in matrices:
        n=len(b);m=cap(b)
        deletions=[cap(np.delete(np.delete(b,i,0),i,1)) for i in range(n)]
        ext,row,count=extension(b)
        c=dict(source=label,source_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
               n=n,M=m,known_exact_M=KNOWN.get(n),deletion_caps=deletions,
               optimal_predecessor_vertices=[i for i,v in enumerate(deletions)
                                             if v==KNOWN.get(n-1)],
               best_extension=ext,E=ext-m,best_row=row,best_row_count=count)
        cases.append(c)
        print(json.dumps(c),flush=True)
out=Path('computations/decisive_independent_nesting_audit_2026_09_06.json')
out.write_text(json.dumps(dict(status='exact finite upper/profile replay; stored lower certificates not rerun',cases=cases),indent=2)+'\n')
