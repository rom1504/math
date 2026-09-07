"""Bounded CP-SAT cap-28 test with two frozen actual order-eight children.

All projective spin-pair inequalities are imposed. UNKNOWN is not an
infeasibility certificate. Returned witnesses receive independent cube checks.
"""
import json
from pathlib import Path
import numpy as np
from ortools.sat.python import cp_model
from flatify_adversary_2026_09_07_hadamard_bridge_phase_search import spins


def run():
    base=Path(__file__).resolve().parent/'results'
    prior=json.loads((base/'flatify_adversary_2026_09_07_frozen_hadamard_trades.json').read_text())
    parent=np.array(prior['parent'],dtype=np.int16)
    a=parent[:8,:8]; d=parent[8:,8:]; x=spins(8)
    ha=np.sum((x@a)*x,axis=1)//2; hd=np.sum((x@d)*x,axis=1)//2
    model=cp_model.CpModel(); z=[model.NewBoolVar('z%d'%i) for i in range(64)]
    old=((1-parent[:8,8:].ravel())//2).tolist()
    for q,value in enumerate(old): model.AddHint(z[q],int(value))
    model.Add(z[0]==int(old[0]))
    for i in range(len(x)):
        for j in range(len(x)):
            f=np.outer(x[i],x[j]).ravel(); budget=28-abs(int(ha[i]+hd[j]))
            linear=sum(int(f[q])*z[q] for q in range(64))
            total=int(f.sum())
            model.Add(linear >= (total-budget+1)//2)
            model.Add(linear <= (total+budget)//2)
    solver=cp_model.CpSolver(); solver.parameters.max_time_in_seconds=120
    solver.parameters.num_search_workers=2
    solver.parameters.random_seed=9070754
    status=solver.Solve(model)
    out=dict(status=solver.StatusName(status),target_cap=28,
             branches=solver.NumBranches(),conflicts=solver.NumConflicts(),
             wall_time=solver.WallTime(),feasible_witness=False)
    if status in (cp_model.OPTIMAL,cp_model.FEASIBLE):
        b=np.array([1-2*solver.Value(v) for v in z],dtype=np.int16).reshape(8,8)
        c=np.block([[a,b],[b.T,d]]); y=spins(16)
        e=np.sum((y@c)*y,axis=1)//2; cap=int(np.max(np.abs(e)))
        assert cap<=28
        out.update(feasible_witness=True,cap=cap,bridge=b.tolist(),parent=c.tolist(),
                   bridge_gram=(b@b.T).tolist())
    (base/'flatify_adversary_2026_09_07_frozen_bridge_cpsat.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True)


if __name__=='__main__': run()
