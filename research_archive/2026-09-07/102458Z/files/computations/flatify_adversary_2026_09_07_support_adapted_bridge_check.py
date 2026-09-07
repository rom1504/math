"""Exact realization and finite angular checks for the all-paired operation."""
import itertools
import json
from fractions import Fraction as Q
from pathlib import Path
import numpy as np
from scipy.linalg import hadamard

records=[]
for k in [4,8]:
    m=k//2; n=m*k
    H=hadamard(k).astype(np.int64); H2=hadamard(2).astype(np.int64)
    active=list(range(0,k,2)); inactive=list(range(1,k,2))
    def frame(port):
        cols=[]
        for j in range(m):
            cols.extend([active[(j-port)%m],inactive[j]])
        return H[:,cols]
    left=[frame(i) for i in range(m)]
    right=[frame((j+1)%m) for j in range(m)]
    C=np.block([[left[i][:,2*j:2*j+2]@H2@right[j][:,2*i:2*i+2].T//2
                 for j in range(m)] for i in range(m)])
    assert np.all(abs(C)==1) and np.array_equal(C@C.T,n*np.eye(n,dtype=np.int64))
    assert int(C.sum())==0
    E=np.repeat(np.eye(n//2,dtype=np.int64),2,axis=0)
    twice=E.T@C@E
    assert np.all(twice%2==0)
    restricted=twice//2
    assert np.array_equal(restricted@restricted.T,(n//2)*np.eye(n//2,dtype=np.int64))
    if k==4:
        values=list(itertools.product([-1,1],repeat=n//2))
        tested=0; maximum=0
        for u in values:
            x=E@np.array(u,dtype=np.int64)
            for v in values:
                y=E@np.array(v,dtype=np.int64)
                if abs(int(x.sum()))!=abs(int(y.sum())): continue
                alpha=Q(abs(int(x.sum())),n); square=alpha*alpha
                bound_squared=Q(n**3,2) if square<=Q(1,2) else 2*square*(1-square)*n**3
                energy=int(x@C@y)
                assert energy*energy<=bound_squared
                tested+=1; maximum=max(maximum,abs(energy))
    else:
        tested=0; maximum=None
    records.append({'k':k,'n':n,'center_bridge':0,'bridge':C.tolist(),
                    'restricted_operator_square':n//2,
                    'exhaustive_equal_overlap_pairs':tested,
                    'restricted_maximum_when_enumerated':maximum})
result={'status':'PASS: EXACT SIGN BRIDGES, CENTER ZERO, RESTRICTED ORTHOGONALITY',
        'records':records,'critical_seed':'sqrt(2+sqrt(2))/4'}
Path('computations/results/flatify_adversary_2026_09_07_support_adapted_bridge_check.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'status':result['status'],'records':[{k:v for k,v in r.items() if k!='bridge'} for r in records]},indent=2))
