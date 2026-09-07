"""Exact full-cut CP-SAT search for an order-11 polarity counterexample.

Any witness is integer-enumerated. UNSAT is a solver report, NOT a portable
proof certificate. The imported M_11=17 lower proof is a separate dependency.
"""
import argparse
import hashlib
import itertools
import json
from pathlib import Path
import time
from ortools.sat.python import cp_model

def run(degree, seconds, seed, log_path):
    n=11; cap=17; reverse_cap=11
    edges=list(itertools.combinations(range(n),2))
    model=cp_model.CpModel()
    neg=[model.NewBoolVar('e_%d_%d'%e) for e in edges]
    model.Add(sum(neg)==19)
    deg=[sum(neg[k] for k,e in enumerate(edges) if i in e) for i in range(n)]
    model.Add(deg[0]==degree)
    for d in deg[1:]: model.Add(d<=degree)
    for j in range(1,n): model.Add(neg[edges.index((0,j))]==int(j<=degree))
    # Remaining permutations may sort degrees separately within each class.
    for group in [list(range(1,degree+1)),list(range(degree+1,n))]:
        for i,j in zip(group,group[1:]): model.Add(deg[i]>=deg[j])
    for size in range(1,6):
        total=size*(n-size)
        lo=max(0,(total-(cap+reverse_cap)//2+1)//2)
        hi=total//2
        for subset in itertools.combinations(range(n),size):
            S=set(subset)
            cut=sum(neg[k] for k,(i,j) in enumerate(edges) if (i in S)!=(j in S))
            model.Add(cut>=lo); model.Add(cut<=hi)
    solver=cp_model.CpSolver()
    solver.parameters.max_time_in_seconds=seconds
    solver.parameters.num_search_workers=1
    solver.parameters.random_seed=seed
    solver.parameters.log_search_progress=True
    solver.parameters.log_to_stdout=False
    start=time.monotonic()
    with log_path.open('w') as stream:
        solver.log_callback=lambda line: (stream.write(line),stream.flush())
        status=solver.Solve(model)
    out=dict(n=n, positive_cap=cap, reverse_cap_upper=reverse_cap,
             fixed_max_negative_degree=degree, seconds=time.monotonic()-start,
             seed=seed, workers=1, status=solver.StatusName(status),
             response_stats=solver.ResponseStats(),
             scope='witness discovery; UNSAT is not independently certified')
    if status in (cp_model.OPTIMAL,cp_model.FEASIBLE):
        values=[1-2*solver.Value(v) for v in neg]
        energies=[]
        for mask in range(1<<(n-1)):
            x=[1]+[1-2*((mask>>i)&1) for i in range(n-1)]
            energies.append(sum(a*x[i]*x[j] for a,(i,j) in zip(values,edges)))
        assert max(energies)==cap and min(energies)>=-reverse_cap
        out.update(edge_signs=values,minimum=min(energies),maximum=max(energies),
                   sha256=hashlib.sha256(bytes(a+1 for a in values)).hexdigest(),
                   integer_verified_states=len(energies))
    return out

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--seconds',type=float,default=1200)
    parser.add_argument('--seed',type=int,default=20260907)
    parser.add_argument('--output',required=True)
    args=parser.parse_args()
    path=Path(args.output); path.parent.mkdir(parents=True,exist_ok=True)
    data=dict(question='Does an order-11 cap-17 signing have polarity gap >=6?',
              imported_minimum_scope='M11=17 from existing separate proof',cases=[])
    for degree in (4,5):
        record=run(degree,args.seconds,args.seed,path.with_suffix('.degree%d.log'%degree))
        data['cases'].append(record)
        path.write_text(json.dumps(data,indent=2)+'\n')
        print(json.dumps(record),flush=True)
