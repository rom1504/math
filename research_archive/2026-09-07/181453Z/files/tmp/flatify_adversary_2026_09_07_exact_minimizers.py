"""Exact switching-reduced full-sign census. No heuristic minimizers.

Gauge first-row signs to +1; x_0=+1 removes global spin symmetry.
Enumerate every remaining signing and every spin, using integer energies.
"""
import json
from pathlib import Path
import numpy as np

OUT = Path(__file__).with_suffix('.json')
results = []
for n in range(3, 8):
    spins = 1 - 2 * ((np.arange(1 << (n-1))[:, None] >> np.arange(n-1)) & 1)
    spins = np.column_stack([np.ones(len(spins), dtype=int), spins])
    pairs = [(i,j) for i in range(1,n) for j in range(i+1,n)]
    tests = np.array([spins[:,i]*spins[:,j] for i,j in pairs], dtype=np.int16)
    star = spins[:,1:].sum(axis=1)
    hist = {}
    examples = {}
    for lo in range(0, 1 << len(pairs), 4096):
        ids = np.arange(lo, min(lo+4096,1 << len(pairs)), dtype=np.int64)
        signs = (1-2*((ids[:,None] >> np.arange(len(pairs))) & 1)).astype(np.int16)
        energies = signs @ tests + star
        p, r = energies.max(axis=1), -energies.min(axis=1)
        for idx,pp,rr in zip(ids,p,r):
            key = (int(pp),int(rr))
            hist[key] = hist.get(key,0)+1
            examples.setdefault(key,int(idx))
    cap = min(max(p,r) for p,r in hist)
    width2 = min(p+r for p,r in hist)
    cap_pairs = sorted(k for k in hist if max(k)==cap)
    width_pairs = sorted(k for k in hist if sum(k)==width2)
    row = dict(n=n,signings=sum(hist.values()),cap=cap,width=width2/2,
               cap_minimizer_pairs=[dict(P=p,R=r,count=hist[(p,r)],id=examples[(p,r)]) for p,r in cap_pairs],
               width_minimizer_pairs=[dict(P=p,R=r,count=hist[(p,r)],id=examples[(p,r)]) for p,r in width_pairs],
               histogram={f'{p},{r}':v for (p,r),v in sorted(hist.items())})
    results.append(row)
    print(json.dumps({k:v for k,v in row.items() if k!='histogram'}),flush=True)
OUT.write_text(json.dumps(results,indent=2)+'\n')
