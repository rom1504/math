"""Exact Gaussian integration; numerical optimization is diagnostic only."""
import itertools
import json
from pathlib import Path
import numpy as np
from scipy.optimize import differential_evolution

t=2*np.sqrt(15.); lam=7.5
R=np.array([[1.,1.],[1.,-1.]])/np.sqrt(2)
T=np.block([[np.zeros((2,2)),R],[R,np.zeros((2,2))]])
A=lam*np.eye(4)-t*T
signs=np.array(list(itertools.product([-1.,1.], repeat=4)))

def kernel(a):
    a=np.asarray(a); V=np.diag(1-a*a)
    P=np.eye(4)+2*V@A
    B=A@np.linalg.inv(P)
    means=signs*a
    return float(np.mean(np.exp(-np.einsum('si,ij,sj->s',means,B,means)))/np.sqrt(np.linalg.det(P)))

corners=[dict(a=list(a),kernel=kernel(a)) for a in itertools.product([0.,1.],repeat=4)]
opt=differential_evolution(lambda a:-kernel(a),[(0.,1.)]*4,seed=41,tol=1e-10,popsize=25,maxiter=500)
result=dict(status='exact integration, uncertified numerical maximum',corners=corners,
            maximum=dict(a=opt.x.tolist(),kernel=-float(opt.fun)),
            semicube_lower=1/16*(1+np.exp(-45+8*np.sqrt(30.))))
print(json.dumps(result,indent=2))
Path('computations/results/flatify_construct_2026_09_07_four_label_kernel.json').write_text(json.dumps(result,indent=2)+'\n')
