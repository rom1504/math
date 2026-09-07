"""Exact pentagon ground-kernel falsifier and interval uniform face gap."""
import itertools
import numpy as np
import sympy as sp
import mpmath as mp
import json
from pathlib import Path

a=np.ones((5,5),dtype=np.int64)-np.eye(5,dtype=np.int64)
for i in range(5):
    a[i,(i+1)%5]=a[(i+1)%5,i]=-1
spins=np.array(list(itertools.product([-1,1],repeat=5)),dtype=np.int64)
energy=np.einsum('bi,ij,bj->b',spins,a,spins)//2
assert set(energy)=={-4,0,4}
assert sum(energy**2)==32*10
ground=spins[energy==4]
triple=np.array([[-1,1,-1,1,1],[-1,1,1,-1,1],[-1,1,-1,-1,1]])
assert np.all(np.einsum('bi,ij,bj->b',triple,a,triple)//2==4)
assert (triple@a@triple.T).tolist()==[[8,4,8],[4,8,8],[8,8,8]]

# Enumerate all positive-ground edges via twice their midpoint.
edges={tuple(x+y) for i,x in enumerate(ground) for y in ground[i+1:] if np.sum(x!=y)==1}
mid=np.array([-1,1,-1,0,1])
generated={tuple(sign*2*np.roll(mid,j)) for sign in [-1,1] for j in range(5)}
assert len(edges)==10 and edges==generated

s=sp.Symbol('s',real=True)
z=sp.Matrix([-1,1,-1,s,1]); aa=sp.Matrix(a)
gram=[]
for j in range(5):
    zj=sp.Matrix([z[(i-j)%5] for i in range(5)])
    gram.append(sp.expand((z.T*aa*zj)[0]))
assert gram==[8,-7-s**2,3+s**2,3+s**2,-7-s**2]

mp.iv.dps=50
iv=mp.iv
phi=(1+iv.sqrt(5))/2
cosh=lambda x:(iv.exp(x)+iv.exp(-x))/2
lam=cosh(iv.mpf(1))-phi*cosh(iv.mpf(7)/8)+cosh(iv.mpf(3)/8)/phi
assert lam < -iv.mpf(728)/10000
print('EXACT_GROUND_KERNEL_AND_UNIFORM_FACE_CERTIFICATE_PASS')
print('lambda_max_over_ground_edges:',lam)

for n,source,key,indices,weights,expected in [
    (10,'exact_m10.json','matrix',[7,8,17,19],[1,-1,-1,1],
     [[26,22,22,2],[22,26,2,-10],[22,2,26,22],[2,-10,22,26]]),
    (14,'conference_completion_m13.json','conference_matrix',[0,1,6],[1,-2,1],
     [[42,38,18],[38,42,38],[18,38,42]])]:
    aa=np.array(json.loads((Path('computations/results')/source).read_text())[key],dtype=np.int64)
    xx=1-2*((np.arange(1<<(n-1))[:,None]>>np.arange(n))&1)
    ee=np.einsum('bi,ij,bj->b',xx,aa,xx)//2
    pp=int(ee.max()); gg=xx[ee==pp][indices]
    gram=gg@aa@gg.T
    assert gram.tolist()==expected
    value=sum(weights[i]*weights[j]*cosh(iv.mpf(int(gram[i,j]))/(2*pp))
              for i in range(len(weights)) for j in range(len(weights)))
    assert value < -iv.mpf(1)/100
    print(json.dumps({'n':n,'P':pp,'ground_spins':gg.tolist(),'Gram':expected,'weights':weights}))
    print('certified negative quadratic:',value)
