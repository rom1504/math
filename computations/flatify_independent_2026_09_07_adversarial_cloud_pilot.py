"""Search rare centers; a heuristic witness is always verified exactly."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import json
from pathlib import Path
import numpy as np
from scipy.linalg import hadamard
rng=np.random.default_rng(202609070921)
records=[]
for k in [8,16,32]:
    m=k//2; n=k*m; f=round(n/10)
    H=hadamard(k).astype(np.int64); H2=hadamard(2).astype(np.int64)
    F=[H[:,rng.permutation(k)]*rng.choice([-1,1],k) for _ in range(m)]
    G=[H[:,rng.permutation(k)]*rng.choice([-1,1],k) for _ in range(m)]
    C=np.block([[F[i][:,2*j:2*j+2]@H2@G[j][:,2*i:2*i+2].T//2 for j in range(m)] for i in range(m)])
    for trial in range(150):
        x=rng.choice([-1,1],n); fields=x@C; y=np.where(fields>=0,1,-1)
        value=int(x@C@y)
        if value/n**1.5<=.74: continue
        X=np.tile(x,(64,1))
        for xx in X: xx[rng.choice(n,f,replace=False)]*=-1
        minimum=10**20; witness=None
        for iteration in range(25):
            fields=-(X@C)*y
            indices=np.argpartition(fields,f-1,axis=1)[:,:f]
            Y=np.tile(y,(64,1)); Y[np.arange(64)[:,None],indices]*=-1
            fields=-(Y@C.T)*x
            indices=np.argpartition(fields,f-1,axis=1)[:,:f]
            X=np.tile(x,(64,1)); X[np.arange(64)[:,None],indices]*=-1
            vals=np.einsum('ij,ij->i',X@C,Y); ix=int(np.argmin(vals))
            if int(vals[ix])<minimum:
                minimum=int(vals[ix]); witness=(X[ix].copy(),Y[ix].copy())
        rec={'n':n,'flips':f,'trial':trial,'target':value,'minimum_center':minimum,'normalized_target':value/n**1.5,'normalized_center':minimum/n**1.5}
        records.append(rec)
        if minimum<=0:
            rec.update({'C':C.tolist(),'F':[a.tolist() for a in F],'G':[a.tolist() for a in G],'x':x.tolist(),'y':y.tolist(),'x0':witness[0].tolist(),'y0':witness[1].tolist()})
            print(json.dumps({key:val for key,val in rec.items() if key not in ['C','F','G','x','y','x0','y0']}),flush=True)
            break
    print('completed order',n,'minimum',min(r['normalized_center'] for r in records if r['n']==n),flush=True)
Path('computations/results/flatify_independent_2026_09_07_adversarial_cloud_pilot.json').write_text(json.dumps(records,indent=2)+'\n')
