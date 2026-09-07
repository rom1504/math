"""Exact finite CP-SAT orientation search at an actual width minimizer."""
import argparse
import itertools
import json
import math
import numpy as np
from ortools.sat.python import cp_model
from decisive_bridge_actual_cavity_psd_checks_2026_09_07 import objects

def run(n,seconds,brute=False):
    edges,chars,signs=objects(n)
    energies=signs@chars.T
    p=energies.max(1)
    r=(-energies).max(1)
    width=(p+r)/2
    ix=int(np.argmin(width))
    a=signs[ix].astype(np.int64)
    model=cp_model.CpModel()
    orient=[model.new_bool_var('e'+str(e)) for e in range(len(edges))]
    bound=model.new_int_var(0,len(edges),'bound')
    constraints=[]
    for subset in range(1,(1<<n)-1):
        S=[i for i in range(n) if subset&(1<<i)]
        T=[i for i in range(n) if not subset&(1<<i)]
        free=[i for i in range(n) if i not in [S[0],T[0]]]
        for bits in itertools.product([-1,1],repeat=len(free)):
            spin=np.ones(n,dtype=np.int64)
            spin[free]=bits
            coefficient=np.zeros(len(edges),dtype=np.int64)
            constant=0
            for e,(i,j) in enumerate(edges):
                v=int(a[e]*spin[i]*spin[j])
                if i in S and j in T:
                    coefficient[e]=v
                elif i in T and j in S:
                    coefficient[e]=-v
                    constant+=v
            expr=constant+sum(int(v)*orient[e] for e,v in enumerate(coefficient) if v)
            model.add(expr<=bound)
            model.add(expr>=-bound)
            constraints.append((constant,coefficient))
    model.minimize(bound)
    solver=cp_model.CpSolver()
    solver.parameters.max_time_in_seconds=seconds
    solver.parameters.num_search_workers=8
    status=solver.solve(model)
    result=dict(n=n,width=float(width[ix]),positive=float(p[ix]),negative=float(r[ix]),status=solver.status_name(status),lower_bound=solver.best_objective_bound,wall_time=solver.wall_time,orientation_target=float(width[ix]/math.sqrt(2)))
    if status in [cp_model.OPTIMAL,cp_model.FEASIBLE]:
        o=np.array([solver.value(v) for v in orient],dtype=np.int64)
        value=max(abs(c+int(v@o)) for c,v in constraints)
        assert value==solver.value(bound)
        result.update(directed_cap=value,lift_cap=4*value,ratio_to_width=value/width[ix],seed_signs=a.tolist(),orientation=o.tolist())
        amat=np.zeros((n,n),dtype=np.int64)
        cmat=np.zeros((n,n),dtype=np.int64)
        for e,(i,j) in enumerate(edges):
            amat[i,j]=amat[j,i]=a[e]
            cmat[i,j]=a[e]*(1-2*o[e])
            cmat[j,i]=-cmat[i,j]
        lift=np.block([[amat,cmat],[-cmat,-amat]])
        permutation=list(range(n,2*n))+list(range(n))
        assert np.array_equal(lift[np.ix_(permutation,permutation)],-lift)
        lifted_spins=np.array([(1,)+s for s in itertools.product([-1,1],repeat=2*n-1)],dtype=np.int64)
        twice_energy=np.sum(lifted_spins*(lifted_spins@lift),axis=1)
        assert int(twice_energy.max())==int((-twice_energy).max())==8*value
        result['lift_identity_full_spin_check']=True
    if brute:
        assert n<=6
        coeff=np.stack([v for c,v in constraints]).astype(np.int16)
        const=np.array([c for c,v in constraints],dtype=np.int16)
        best=len(edges)
        count=1<<len(edges)
        shifts=np.arange(len(edges),dtype=np.int64)
        for start in range(0,count,512):
            words=np.arange(start,min(start+512,count),dtype=np.int64)
            orientations=((words[:,None]>>shifts)&1).astype(np.int16)
            values=np.max(np.abs(orientations@coeff.T+const),axis=1)
            best=min(best,int(values.min()))
        assert best==int(solver.best_objective_bound)==int(solver.objective_value)
        result['independent_exhaustive_orientation_count']=count
        result['independent_exact_minimum']=best
    print(json.dumps(result),flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--n',type=int,default=6)
    parser.add_argument('--seconds',type=float,default=30)
    parser.add_argument('--brute',action='store_true')
    args=parser.parse_args()
    run(args.n,args.seconds,args.brute)
