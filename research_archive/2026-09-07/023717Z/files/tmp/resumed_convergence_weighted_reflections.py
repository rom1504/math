"""Numerical exploration only: two weighted binary averaging reflections."""
import itertools
import numpy as np
from scipy.optimize import differential_evolution

B=np.array([[-1.,2.],[2.,4.]])
profiles=np.array(list(itertools.product([-1.,1.],repeat=8))).reshape(-1,4,2)

def value(z,ret=False):
    p=np.array([z[0],1-z[0]])
    q=np.array([z[1],1-z[1]])
    mu=np.kron(p,q)
    U=np.kron(2*np.tile(p,(2,1))-np.eye(2),2*np.tile(q,(2,1))-np.eye(2))
    fields=np.einsum('ab,tbi,ij->taj',U,profiles,B)
    vals=(np.abs(fields).sum(axis=2)*mu).sum(axis=1)/2
    k=vals.argmax()
    return (vals[k],profiles[k],np.sign(fields[k])) if ret else -vals[k]

res=differential_evolution(value,[(0.001,.5),(.001,.5)],tol=1e-11,popsize=20,
                           maxiter=500,seed=20260906,polish=True)
print(res.x,-res.fun, 'target',5/np.sqrt(2),flush=True)
print(value(res.x,True),flush=True)
