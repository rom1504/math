#!/usr/bin/env python3
"""Find diverse globally optimal small children, without claiming class exhaustiveness."""
import argparse
import itertools
import json
from pathlib import Path
import time

import networkx as nx
import numpy as np
from ortools.sat.python import cp_model

ROOT=Path(__file__).resolve().parents[1]


def rooted_graph(a,root,polarity=1):
    n=len(a);p=[root]+[i for i in range(n) if i!=root]
    b=a[np.ix_(p,p)]
    s=np.r_[1,b[0,1:]]
    b=polarity*s[:,None]*b*s[None,:]
    g=nx.Graph();g.add_nodes_from(range(n-1))
    g.add_edges_from((i-1,j-1) for i in range(1,n) for j in range(i+1,n) if b[i,j]>0)
    return g


def equivalent(a,b):
    ga=rooted_graph(a,0)
    signature=sorted(dict(ga.degree()).values())
    for root in range(len(a)):
        for polarity in (1,-1):
            gb=rooted_graph(b,root,polarity)
            if sorted(dict(gb.degree()).values())==signature and nx.is_isomorphic(ga,gb):return True
    return False


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--order',type=int,required=True)
    parser.add_argument('--cap',type=int,required=True)
    parser.add_argument('--samples',type=int,default=100)
    parser.add_argument('--time-limit',type=float,default=5)
    parser.add_argument('--seed',type=int,default=20260918)
    parser.add_argument('--tag',default='sampled')
    args=parser.parse_args();n=args.order
    rng=np.random.default_rng(args.seed)
    edges=list(itertools.combinations(range(n),2));free=[e for e in edges if e[0]]
    codes=np.arange(1<<(n-1),dtype=np.uint64)
    spins=np.c_[np.ones(len(codes),dtype=np.int16),1-2*((codes[:,None]>>np.arange(n-1,dtype=np.uint64))&1).astype(np.int16)]
    model=cp_model.CpModel();variables=[model.new_bool_var(f'a{i}_{j}') for i,j in free]
    products=np.asarray([spins[:,i]*spins[:,j] for i,j in free]).T
    constant=np.asarray([spins[:,i]*spins[:,j] for i,j in edges]).sum(axis=0)
    for co,c in zip(products,constant):
        expression=sum(int(v)*b for v,b in zip(co,variables))
        model.add(expression>=int((c-args.cap+1)//2));model.add(expression<=int((c+args.cap)//2))
    initial=np.asarray(json.loads((ROOT/f'computations/results/exact_m{n}.json').read_text())['matrix'],dtype=np.int16)
    seeds=[dict(label=f'm{n}_sample0',matrix=initial.tolist(),cap=args.cap,source=f'computations/results/exact_m{n}.json')]
    reps=[initial]
    log=[];started=time.time()
    for sample in range(args.samples):
        coefficients=rng.integers(-100,101,len(variables))
        model.maximize(sum(int(c)*b for c,b in zip(coefficients,variables)))
        solver=cp_model.CpSolver();solver.parameters.max_time_in_seconds=args.time_limit
        solver.parameters.num_search_workers=1;solver.parameters.random_seed=int(rng.integers(1,2**30))
        status=solver.solve(model)
        row=dict(sample=sample,status=solver.status_name(status),elapsed=solver.wall_time)
        if status in (cp_model.OPTIMAL,cp_model.FEASIBLE):
            a=np.ones((n,n),dtype=np.int16);np.fill_diagonal(a,0)
            bits=[solver.value(b) for b in variables]
            for (i,j),b in zip(free,bits):a[i,j]=a[j,i]=1-2*b
            actual=int(np.max(np.abs(np.einsum('bi,ij,bj->b',spins,a,spins)//2)))
            assert actual<=args.cap
            cls=next((i for i,b in enumerate(reps) if equivalent(a,b)),None)
            if cls is None:
                cls=len(reps);reps.append(a)
                seeds.append(dict(label=f'm{n}_sample{cls}',matrix=a.tolist(),cap=actual,
                    source='CP-SAT feasible cap-constrained sample; exact child cap checked; class list not exhaustive',
                    random_seed=args.seed,sample=sample,objective_coefficients=coefficients.tolist()))
                print(f'NEW order={n} class={cls} sample={sample} cap={actual} elapsed={time.time()-started:.3f}',flush=True)
            row.update(class_id=cls,cap=actual)
            # Avoid the exact same root-gauged solution in subsequent samples.
            model.add(sum((1-b if bit else b) for b,bit in zip(variables,bits))>=1)
        log.append(row)
        (ROOT/f'computations/results/twisted_chiral_2026_09_18_{args.tag}_m{n}_seeds.json').write_text(json.dumps(seeds,indent=2)+'\n')
        (ROOT/f'computations/results/twisted_chiral_2026_09_18_{args.tag}_m{n}_sampling.json').write_text(json.dumps(log,indent=2)+'\n')
    print(f'DONE order={n} distinct_classes={len(seeds)} elapsed={time.time()-started:.3f}',flush=True)


if __name__=='__main__':main()
