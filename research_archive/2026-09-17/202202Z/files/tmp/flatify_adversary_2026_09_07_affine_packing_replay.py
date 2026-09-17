"""Exact independent pentagon/affine-spike checks, including arbitrary D_i."""
import itertools
import json
from pathlib import Path
import numpy as np

rng=np.random.default_rng(202609070610)
A=np.ones((5,5),dtype=np.int64)-np.eye(5,dtype=np.int64)
for i in range(5):A[i,(i+1)%5]=A[(i+1)%5,i]=-1
ternary=np.array(list(itertools.product([-1,0,1],repeat=5)),dtype=np.int64)
energy=np.einsum('bi,ij,bj->b',ternary,A,ternary)//2
assert np.all(np.abs(energy)<=np.sum(ternary*ternary,axis=1))
full=ternary[np.all(ternary!=0,axis=1)]
assert np.max(np.abs(np.einsum('bi,ij,bj->b',full,A,full)))==8
patterns={}
for polarity in [-1,1]:
    patterns[polarity]=[]
    for S in itertools.combinations(range(5),3):
        for s in itertools.product([-1,1],repeat=3):
            z=np.zeros(5,dtype=np.int64);z[list(S)]=s
            if polarity*int(z@A@z//2)==3:
                patterns[polarity].append((S,z));break
    assert len(patterns[polarity])==5
    assert all(sum(i in S for S,_ in patterns[polarity])==3 for i in range(5))

def mul4(a,b):
    a0,a1=a&1,a>>1;b0,b1=b&1,b>>1
    return (a0*b0^a1*b1)+2*(a0*b1^a1*b0^a1*b1)

q=4
lines=[]
for slope in [None,0,1,2,3]:
    family=[]
    for c in range(q):
        family.append([(c,y) for y in range(q)] if slope is None else [(x,mul4(slope,x)^c) for x in range(q)])
    lines.append(family)
for i in range(5):
    assert len(set(p for L in lines[i] for p in L))==q*q
    for j in range(i):
        assert all(len(set(L)&set(M))==1 for L in lines[i] for M in lines[j])
# The first attempt used int.bit_count, unavailable in this Python 3.9 venv;
# changed to bin(...).count('1') before any claimed test completion.
H=np.array([[(-1)**bin(a&b).count('1') for b in range(q)] for a in range(q)],dtype=np.int64)
cases=[]
for trial in range(20):
    line_h={}
    V=[]
    for i in range(5):
        v=np.zeros((q*q,q*q),dtype=np.int64)
        for l,L in enumerate(lines[i]):
            h=H[rng.permutation(q)][:,rng.permutation(q)]*rng.choice([-1,1],q)[:,None]
            line_h[i,l]=h
            for c,(x,y) in enumerate(L):v[l*q:(l+1)*q,q*x+y]=h[:,c]
        assert np.array_equal(v@v.T,q*np.eye(q*q,dtype=np.int64))
        V.append(v)
    F=np.zeros((5*q*q,5*q*q),dtype=np.int64)
    for i in range(5):
        d=rng.choice([-1,1],(q*q,q*q));d=np.triu(d,1);d=d+d.T
        F[i*q*q:(i+1)*q*q,i*q*q:(i+1)*q*q]=d
        for j in range(i):
            c=A[i,j]*(V[i]@V[j].T)
            assert np.all(np.abs(c)==1)
            F[i*q*q:(i+1)*q*q,j*q*q:(j+1)*q*q]=c
            F[j*q*q:(j+1)*q*q,i*q*q:(i+1)*q*q]=c.T
    assert np.array_equal(F,F.T) and np.all(np.diag(F)==0)
    assert np.all(np.abs(F+np.eye(len(F),dtype=np.int64))==1)
    for polarity in [-1,1]:
        edges=[];degrees={};codegrees={};total=0
        for p in itertools.product(range(q),repeat=2):
            for S,z in patterns[polarity]:
                vertices=tuple((i,next(l for l,L in enumerate(lines[i]) if p in L)) for i in S)
                cols=[];value=3*q*q
                for i,l in vertices:
                    h=line_h[i,l][:,lines[i][l].index(p)]
                    ix=np.arange(i*q*q+l*q,i*q*q+(l+1)*q)
                    value+=polarity*int(h@F[np.ix_(ix,ix)]@h//2)
                    cols.append((ix,z[i]*h))
                    degrees[i,l]=degrees.get((i,l),0)+1
                for pair in itertools.combinations(vertices,2):codegrees[pair]=codegrees.get(pair,0)+1
                edges.append((vertices,cols,value))
                total+=value
        assert len(edges)==5*q*q and set(degrees.values())=={3*q}
        assert max(codegrees.values())<=3
        assert total==15*q**4
        used=set();selected=[]
        for idx in rng.permutation(len(edges)):
            edge=edges[idx]
            if not used.intersection(edge[0]):selected.append(edge);used.update(edge[0])
        covariance=np.eye(len(F),dtype=np.int64)
        for vertices,cols,value in selected:
            vector=np.zeros(len(F),dtype=np.int64)
            for ix,h in cols:vector[ix]=h;covariance[ix,ix]=0
            covariance+=np.outer(vector,vector)
        twice=int(np.sum(F*covariance))
        assert twice%2==0
        measured=polarity*twice//2
        predicted=sum(edge[2] for edge in selected)
        assert measured==predicted
        cases.append(dict(trial=trial,polarity=polarity,edges=len(edges),degree=3*q,
                          maximum_codegree=max(codegrees.values()),total_weight=total,
                          matching_size=len(selected),mean_signed_energy=measured))
out=dict(status='exact independent PASS',seed=A.tolist(),ternary_states=len(ternary),
         parent_order=5*q*q,cases=cases)
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(dict(status=out['status'],cases=len(cases),ternary_states=len(ternary),parent_order=len(F))))
