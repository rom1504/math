"""Finite numerical check of the convex-orbit difference-norm obstruction.

n=4, two exact optimal order-two children. The orbit is all signed perfect
matchings, scaled by sqrt(3). Enumerates every parent modulo switching.
"""
import itertools
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog

n=4
edges=list(itertools.combinations(range(n),2))
spins=np.array([(1,)+x for x in itertools.product([-1,1],repeat=n-1)])
tests=np.array([spins[:,i]*spins[:,j] for i,j in edges]).T
matchings=[[(0,1),(2,3)],[(0,2),(1,3)],[(0,3),(1,2)]]
vertices=[]
for matching in matchings:
    for signs in itertools.product([-1,1],repeat=2):
        v=np.zeros(len(edges))
        for e,s in zip(matching,signs):v[edges.index(e)]=s*np.sqrt(3)
        vertices.append(v)
vertices=np.array(vertices).T
orbit_energy=tests@vertices
out=[]
for tail in itertools.product([-1,1],repeat=3):
    c=np.array([1,1,1,*tail])
    h=tests@c
    matrix=np.vstack([np.column_stack([-orbit_energy,-np.ones(len(spins))]),
                      np.column_stack([orbit_energy,-np.ones(len(spins))])])
    rhs=np.r_[-h,h]
    objective=np.r_[np.zeros(vertices.shape[1]),1]
    result=linprog(objective,A_ub=matrix,b_ub=rhs,
                   A_eq=np.array([np.r_[np.ones(vertices.shape[1]),0]]),b_eq=[1],
                   bounds=[(0,None)]*len(objective),method='highs')
    assert result.success
    d=vertices@result.x[:-1]
    discrepancy=float(np.max(np.abs(tests@(c-d))))
    bound=(1-1/np.sqrt(3))*4*np.sqrt(3)/(4*np.sqrt(2))
    assert discrepancy+1e-8>=bound
    out.append(dict(parent=c.tolist(),parent_cap=int(np.max(np.abs(h))),
                    best_discrepancy=discrepancy,bound=bound,
                    hull_weights=result.x[:-1].tolist(),D=d.tolist()))
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
