"""Lower-bound pilot: exact-cardinality alternating local cloud maximization."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import json
import numpy as np
from scipy.linalg import hadamard

rng=np.random.default_rng(202609070919)
records=[]
for k in [8,16,32]:
    m=k//2; n=m*k; flips=round(n/10)
    base=hadamard(k).astype(np.int64)
    left=[base[:,rng.permutation(k)]*rng.choice([-1,1],size=k) for _ in range(m)]
    right=[base[:,rng.permutation(k)]*rng.choice([-1,1],size=k) for _ in range(m)]
    C=np.empty((n,n),dtype=np.int64)
    H2=np.array([[1,1],[1,-1]],dtype=np.int64)
    for i in range(m):
        for j in range(m):
            C[i*k:(i+1)*k,j*k:(j+1)*k]=(left[i][:,2*j:2*j+2]@H2@right[j][:,2*i:2*i+2].T)//2
    assert np.all(np.abs(C)==1)
    assert np.array_equal(C@C.T,n*np.eye(n,dtype=np.int64))
    for trial in range(12):
        x0=rng.choice([-1,1],size=n); y0=rng.choice([-1,1],size=n)
        X=np.tile(x0,(128,1))
        for row in X:
            row[rng.choice(n,flips,replace=False)]*=-1
        best=-10**20; bestxy=None
        for iteration in range(40):
            fields=(X@C)*y0
            indices=np.argpartition(fields,flips-1,axis=1)[:,:flips]
            Y=np.tile(y0,(128,1)); Y[np.arange(128)[:,None],indices]*=-1
            fields=(Y@C.T)*x0
            indices=np.argpartition(fields,flips-1,axis=1)[:,:flips]
            X=np.tile(x0,(128,1)); X[np.arange(128)[:,None],indices]*=-1
            values=np.einsum('ij,ij->i',X@C,Y)
            ix=int(np.argmax(values))
            if int(values[ix])>best:
                best=int(values[ix]); bestxy=(X[ix].copy(),Y[ix].copy())
        rec={'k':k,'n':n,'flips':flips,'trial':trial,'center':int(x0@C@y0),
             'best':best,'normalized_center':float(x0@C@y0/n**1.5),'normalized_best':float(best/n**1.5)}
        records.append(rec)
        print(json.dumps(rec),flush=True)
        if abs(rec['normalized_center'])<.03 and rec['normalized_best']>.74:
            rec.update({'C':C.tolist(),'x0':x0.tolist(),'y0':y0.tolist(),'x':bestxy[0].tolist(),'y':bestxy[1].tolist()})
Path('computations/results/flatify_independent_2026_09_07_local_cloud_pilot.json').write_text(json.dumps(records,indent=2)+'\n')
