"""Bounded MILP: exact cap-28 feasibility with frozen actual m8 children.

Unlike the Hadamard orbit search, bridge entries are unrestricted signs.
Solver failure/timeout is not a mathematical infeasibility certificate.
Any returned feasible bridge is independently checked on the full cube.
"""
import json
from pathlib import Path
import numpy as np
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import csc_matrix
from flatify_adversary_2026_09_07_hadamard_bridge_phase_search import spins


def run():
    base=Path(__file__).resolve().parent/'results'
    prior=json.loads((base/'flatify_adversary_2026_09_07_frozen_hadamard_trades.json').read_text())
    seed=np.array(prior['parent'],dtype=np.int16); n=8; target_cap=28
    a=seed[:n,:n]; d=seed[n:,n:]; x=spins(n)
    ha=np.sum((x@a)*x,axis=1)//2; hd=np.sum((x@d)*x,axis=1)//2
    budgets=target_cap-np.abs(ha[:,None]+hd[None,:])
    feature=np.einsum('xi,yj->xyij',x,x).reshape(len(x)**2,n*n).astype(float)
    sums=feature.sum(axis=1)
    lo=(sums-budgets.ravel())/2; hi=(sums+budgets.ravel())/2
    old=(1-seed[:n,n:].ravel())/2
    objective=1-2*old
    lower=np.zeros(n*n); upper=np.ones(n*n)
    # Whole-bridge negation leaves its cap contribution unchanged.
    lower[0]=upper[0]=old[0]
    result=milp(objective,integrality=np.ones(n*n),bounds=Bounds(lower,upper),
        constraints=LinearConstraint(csc_matrix(feature),lo,hi),
        options=dict(time_limit=100.,mip_rel_gap=0.))
    output=dict(status=int(result.status),message=str(result.message),
        target_cap=target_cap,child_caps=[int(np.max(np.abs(ha))),int(np.max(np.abs(hd)))],
        solver_fun=None if result.fun is None else float(result.fun),
        feasible_witness=False)
    if result.x is not None and np.max(np.abs(result.x-np.rint(result.x)))<1e-6:
        b=(1-2*np.rint(result.x).astype(np.int16)).reshape(n,n)
        parent=np.block([[a,b],[b.T,d]])
        xp=spins(16)
        energies=np.sum((xp@parent)*xp,axis=1)//2
        cap=int(np.max(np.abs(energies)))
        if cap<=target_cap:
            output.update(feasible_witness=True,cap=cap,bridge=b.tolist(),parent=parent.tolist(),
                          changed_from_hadamard=int(np.sum(b!=seed[:n,n:])),
                          bridge_gram=(b@b.T).tolist())
    print(json.dumps({k:v for k,v in output.items() if k not in ['bridge','parent','bridge_gram']}),flush=True)
    return output


if __name__=='__main__':
    result=run()
    target=Path(__file__).resolve().parent/'results'/'flatify_adversary_2026_09_07_frozen_bridge_milp.json'
    target.write_text(json.dumps(result,indent=2)+'\n')
