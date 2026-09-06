"""Exact cosquare full-sign screening; SDP numbers are diagnostic only."""
import argparse
import itertools
import json
import math
import numpy as np
import cvxpy as cp


def run(n, limit):
    positions = list(zip(*np.triu_indices(n)))
    spins = np.array(list(itertools.product((-1, 1), repeat=n)), dtype=np.int16)
    groups = {}
    for code in range(1 << (len(positions)-1)):
        a = np.ones((n, n), dtype=np.int16)
        for bit, (i,j) in enumerate(positions[1:]):
            a[i,j] = a[j,i] = 1-2*((code >> bit)&1)
        square = a @ a
        beta = int(np.max(np.sum(np.abs(spins @ a), axis=1)))
        key = square.astype(np.int8).tobytes()
        if key not in groups:
            groups[key] = [beta, code, beta, code]
        else:
            item = groups[key]
            if beta < item[0]: item[0], item[1] = beta, code
            if beta > item[2]: item[2], item[3] = beta, code
    candidates = []
    for key, item in groups.items():
        if item[2]>item[0]:
            candidates.append((item[2]-item[0], key, item))
    candidates.sort(reverse=True)
    d = cp.Variable((n,n), symmetric=True)
    aparam = cp.Parameter((n,n), symmetric=True)
    t = cp.Variable()
    constraints = [d-aparam >> 0, d+aparam >> 0]
    constraints += [cp.sum(cp.multiply(np.outer(x,x),d)) <= 2*t for x in spins]
    problem = cp.Problem(cp.Minimize(t), constraints)
    results = []
    for spread, square, item in candidates[:limit]:
        low_beta, low_code, best_beta, best_code = item
        a = np.ones((n,n),dtype=np.int16)
        for bit,(i,j) in enumerate(positions[1:]):
            a[i,j]=a[j,i]=1-2*((low_code>>bit)&1)
        aparam.value=a
        tv=problem.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_gap_rel=1e-9,tol_feas=1e-9,max_iter=200)
        results.append(dict(square=list(map(int,np.frombuffer(square,dtype=np.int8))),
                            low_code=low_code,high_code=best_code,
                            low_beta=low_beta,high_beta=best_beta,
                            T_low_numeric=tv,separation_numeric=best_beta/2-tv))
    print(json.dumps(dict(n=n,full_signs_enumerated=1<<(len(positions)-1),
                          square_groups=len(groups),beta_varying_groups=len(candidates),
                          checked=results),sort_keys=True))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--n',type=int,default=5)
    p.add_argument('--limit',type=int,default=12)
    args=p.parse_args()
    run(args.n,args.limit)
