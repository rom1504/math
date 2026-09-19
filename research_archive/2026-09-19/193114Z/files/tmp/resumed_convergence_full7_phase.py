"""Exhaustive inputs for full seven-point phase family; exploratory floats."""
import numpy as np
from scipy.optimize import minimize_scalar
from resumed_convergence_gauss_phase_search import phase_operator

v = 7
bits = np.arange(2**(2*v-1), dtype=np.uint32)[:,None]
packed = 1 - 2*((bits >> np.arange(2*v)) & 1).astype(np.int8)
F = packed.reshape((-1,v,2)).astype(float)
B = np.array([[-1.,2.],[2.,4.]])


def best(theta):
    U = phase_operator(1,theta)+np.ones((v,v))/v
    values = np.abs(np.einsum('ab,nbc,cd->nad',U,F,B,optimize=True)).sum(axis=(1,2))/(2*v)
    idx = np.argmax(values)
    return values[idx], idx


record = (-1,None,None)
for theta in np.linspace(0,np.pi,501):
    val, idx = best(theta)
    if val>record[0]+1e-10:
        record=val,idx,theta
        print('record', record, flush=True)
f=F[record[1]]
fit = minimize_scalar(lambda th: -np.abs((phase_operator(1,th)+np.ones((v,v))/v)@f@B).sum()/(2*v),
                      bounds=(record[2]-.05,record[2]+.05),method='bounded')
print('BEST',-fit.fun,fit.x,flush=True)
print(f.astype(int).tolist(),flush=True)
