"""All Boolean profiles on full Paley three-cell rotation times Q_p."""
import numpy as np
from scipy.optimize import differential_evolution

weights3=np.array([1.,3.,3.])/7
A3=np.array([[0.,-3.,3.],[1.,0.,-1.],[-1.,1.,0.]])/np.sqrt(7)
E3=np.tile(weights3,(3,1))
assert np.max(np.abs(A3@A3+np.eye(3)-E3))<1e-14
num=6
codes=np.arange(2**(2*num-1),dtype=np.uint32)[:,None]
F=(1-2*((codes>>np.arange(2*num))&1).astype(np.int8)).reshape((-1,num,2)).astype(float)
B=np.array([[-1.,2.],[2.,4.]])
FB=F@B
record=(-1,None,None)


def objective(parameters, verbose=False):
    global record
    theta,p=parameters
    U3=E3+np.cos(theta)*(np.eye(3)-E3)+np.sin(theta)*A3
    Q2=2*np.tile([p,1-p],(2,1))-np.eye(2)
    U=np.kron(Q2,U3)
    weights=np.kron([p,1-p],weights3)
    fields=np.einsum('ab,nbc->nac',U,FB,optimize=True)
    values=.5*np.einsum('nac,a->n',np.abs(fields),weights)
    i=np.argmax(values)
    if values[i]>record[0]+1e-9:
        record=(values[i],parameters.copy(),i)
        print('record',record,flush=True)
    return -values[i]


fit=differential_evolution(objective,[(0,2*np.pi),(0,1)],popsize=20,maxiter=150,
                           seed=712,polish=True,tol=1e-11)
print('BEST',record,flush=True)
print(F[record[2]].astype(int).tolist(),flush=True)
