"""Independent exact quadratic MUB witness/orbit reconstruction.

Uses only the candidate basis generator, not the proposing witness code.
Checks relative Hadamards really have quadratic column phases and all rows
are distinct signed characters times that phase; uses u=3 and largest v.
"""
import itertools
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'computations'))
from flatify_construct_2026_09_07_kerdock_diagnostic import bases16

B=bases16();t=16
A=np.array([[0,1,1,1,1],[1,0,1,1,-1],[1,1,0,-1,-1],[1,1,-1,0,1],[1,-1,-1,1,0]],dtype=np.int64)
patterns=np.array([[-1,1,1,-1],[-1,1,-1,-1],[-1,1,-1,-1],[0,0,2,0],[0,0,2,0]],dtype=np.int64)
assert sum(int(A[i,j]*(patterns[i]@patterns[j])) for i in range(5) for j in range(i+1,5))==20
parity=lambda n:bin(int(n)).count('1')%2

def check_quadratic(signs):
    f=(signs<0).astype(int)
    linear=[f[1<<i]^f[0] for i in range(4)]
    polar={(i,j):f[(1<<i)|(1<<j)]^f[0]^linear[i]^linear[j] for i in range(4) for j in range(i+1,4)}
    for x in range(t):
        v=int(f[0])
        for i in range(4):
            if (x>>i)&1:v^=int(linear[i])
        for (i,j),a in polar.items():
            if (x>>i)&1 and (x>>j)&1:v^=int(a)
        assert v==f[x]
    return f

rng=np.random.default_rng(202609070625)
rows=[];orbit_rows=[]
orbit_cases={(0,1,2,3,4),(1,2,3,4,5),(0,2,4,6,8)}
for ids in itertools.combinations(range(9),5):
    ref=next(i for i in range(9) if i not in ids)
    relative=[];phases=[]
    for i in ids:
        product=B[i]@B[ref].T
        assert np.all(np.abs(product)==4)
        R=product//4
        phase=R[0]
        phases.append(check_quadratic(phase));relative.append(R)
        labels=[]
        for row in R:
            ratio=row*phase
            normalized=ratio*ratio[0]
            label=sum((normalized[1<<a]<0)*(1<<a) for a in range(4))
            assert all(normalized[x]==(-1)**parity(label&x) for x in range(t))
            labels.append(label)
        assert sorted(labels)==list(range(t))
    form=lambda i,u,v:int(phases[i][0]^phases[i][u]^phases[i][v]^phases[i][u^v])
    u=3
    choices=[v for v in range(t) if (form(0,u,v)^form(1,u,v))==1 and (form(0,u,v)^form(2,u,v))==1]
    assert len(choices)==t//4
    v=max(choices);c=form(1,u,v)
    w=np.zeros((5,t),dtype=np.int64)
    adjusted=patterns.copy();adjusted[:,3]*=(-1)**c
    w[:,[0,u,v,u^v]]=adjusted
    product=np.array([R@z for R,z in zip(relative,w)])
    assert np.all(product%2==0)
    x=product//2
    assert np.all(np.abs(x)==1)
    C=np.block([[np.zeros((t,t),dtype=np.int64) if i==j else A[i,j]*(B[ids[i]]@B[ids[j]].T//4)
                 for j in range(5)] for i in range(5)])
    assert int(x.ravel()@C@x.ravel()//2)==320
    rows.append(dict(ids=ids,reference=ref,u=u,v=v,phase_bit=c,witness=x.ravel().tolist()))
    if ids in orbit_cases:
        samples=[]
        for p in range(t):
            for ell in range(t):
                tw=np.array([[(-1)**parity(ell&z)*w[i,z^p] for z in range(t)] for i in range(5)],dtype=np.int64)
                assert np.array_equal(tw@tw.T,w@w.T)
                product=np.array([R@z for R,z in zip(relative,tw)])
                assert np.all(product%2==0)
                xx=product//2
                assert np.all(np.abs(xx)==1)
                assert int(xx.ravel()@C@xx.ravel()//2)==320
                samples.append(xx)
        samples=np.array(samples)
        for i in range(5):assert np.array_equal(samples[:,i].T@samples[:,i],256*np.eye(t,dtype=np.int64))
        completions=[]
        for _ in range(10):
            F=C.copy()
            for i in range(5):
                D=np.triu(rng.choice([-1,1],(t,t)),1);D=D+D.T
                F[i*t:(i+1)*t,i*t:(i+1)*t]=D
            xs=samples.reshape(256,5*t)
            twice=np.einsum('bi,ij,bj->b',xs,F,xs)
            assert np.all(twice%2==0)
            energies=twice//2
            assert int(energies.sum())==256*320
            assert energies.max()>=320
            completions.append(dict(minimum=int(energies.min()),maximum=int(energies.max()),sum=int(energies.sum())))
        orbit_rows.append(dict(ids=ids,orbit_size=256,covariance='256I exact',completions=completions))
out=dict(status='independent exact PASS',all_subsets=len(rows),witness_energy=320,
         signed_character_factorizations=5*len(rows),rows=rows,orbit_cases=orbit_rows)
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(dict(status=out['status'],subsets=len(rows),orbit_cases=len(orbit_rows),arbitrary_completions=30)))
