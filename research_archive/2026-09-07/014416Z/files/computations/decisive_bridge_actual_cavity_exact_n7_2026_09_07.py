"""Exact global optimizer and conditional-PSD counterexample at lambda=log(2)."""
import itertools
from fractions import Fraction
import numpy as np

n=7
edges=list(itertools.combinations(range(n),2))
spins=np.array([(1,)+x for x in itertools.product((-1,1),repeat=n-1)],dtype=np.int64)
chars=np.array([spins[:,i]*spins[:,j] for i,j in edges]).T
free=[q for q,(i,j) in enumerate(edges) if i!=0]
signs=np.ones((2**len(free),len(edges)),dtype=np.int64)
signs[:,free]=np.array(list(itertools.product((-1,1),repeat=len(free))))
energy=signs@chars.T
ab=np.abs(energy)
# Twice the cosh sum times 2^21, entirely integral and below int64 overflow.
weight=np.left_shift(np.int64(1),21+ab)+np.left_shift(np.int64(1),21-ab)
Z=weight.sum(axis=1)
minimum=int(Z.min())
minimizers=np.flatnonzero(Z==minimum)
print('Exact minimizer count:',len(minimizers),'minimum scaled partition:',minimum)

best=None
for ix in minimizers:
    E=energy[ix]
    cs=[]
    for e in range(len(edges)):
        cav=E-signs[ix,e]*chars[:,e]
        plus=np.left_shift(np.int64(1),20+cav)
        minus=np.left_shift(np.int64(1),20-cav)
        den=int((plus+minus).sum())
        num=int(((plus-minus)*chars[:,e]).sum())
        assert int(signs[ix,e])*num<=0
        cs.append(Fraction(abs(num),den))
    K=np.zeros((n,n))
    for c,(i,j) in zip(cs,edges):
        K[i,j]=K[j,i]=float(c)
    for a in range(1,n//2+1):
        for block in itertools.combinations(range(n),a):
            z=np.full(n,-a,dtype=np.int64)
            z[list(block)]=n-a
            form=float(z@K@z+z@z)
            if best is None or form<best[0]:
                best=(form,int(ix),block,z.copy(),cs)

_,ix,block,z,cs=best
form=sum(Fraction(2*int(z[i])*int(z[j]))*c for c,(i,j) in zip(cs,edges))
diag=int(z@z)
assert form+diag<0
print('Witness switching-class index:',ix)
print('Upper-triangle signs:',signs[ix].tolist())
print('Zero-sum integer block vector:',z.tolist())
print('Cavity absolute correlations:',[str(c) for c in cs])
print('Exact z^T K z:',form)
print('Exact z^T (I+K) z:',form+diag,'decimal:',float(form+diag))
print('PASS: actual globally minimal partition, but diagonal-repaired conditional PSD fails.')
