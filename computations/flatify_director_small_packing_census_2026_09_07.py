"""Exact rational support-packing optima on the graph-atlas switching cover.

Enumeration completeness through order 8 uses NetworkX's documented atlas of
all unlabeled graphs through seven vertices. Every displayed cap and every
packing primal/dual certificate is independently checked over integers or
Fractions. LP floating point proposes certificates; it does not validate them.
This is finite evidence, not an asymptotic sign-matrix theorem.
"""
import argparse
from fractions import Fraction as F
import hashlib
import itertools
import json
from pathlib import Path
import time

import networkx as nx
import numpy as np
from scipy.optimize import linprog


def packing(values, incidence):
    n = incidence.shape[0]
    result = linprog(-values.astype(float), A_ub=incidence,
                     b_ub=np.ones(n), bounds=(0, None), method="highs")
    assert result.success
    primal = [(j, F(float(v)).limit_denominator(1000000))
              for j, v in enumerate(result.x) if v > 1e-9]
    dual = [F(float(v)).limit_denominator(1000000)
            for v in -result.ineqlin.marginals]
    assert all(v >= 0 for _, v in primal) and all(v >= 0 for v in dual)
    assert all(sum(v*int(incidence[i,j]) for j,v in primal) <= 1
               for i in range(n))
    assert all(sum(dual[i]*int(incidence[i,j]) for i in range(n)) >= int(values[j])
               for j in range(len(values)))
    lower = sum(v*int(values[j]) for j,v in primal)
    upper = sum(dual)
    assert lower == upper, (lower, upper)
    return dict(value=str(upper),
                primal=[dict(support=j+1, weight=str(v), energy=int(values[j]))
                        for j,v in primal], dual=list(map(str,dual)))


def run(max_n):
    atlas = nx.graph_atlas_g()
    counts = [sum(len(g)==k for g in atlas) for k in range(8)]
    assert counts == [1,1,2,4,11,34,156,1044]
    records=[]
    for n in range(3,max_n+1):
        start=time.monotonic()
        z=np.array(list(itertools.product([-1,0,1],repeat=n)),dtype=np.int16)
        supports=((z!=0)*(1<<np.arange(n))).sum(axis=1)
        pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
        products=np.array([z[:,i]*z[:,j] for i,j in pairs],dtype=np.int16)
        full=supports==(1<<n)-1
        incidence=((np.arange(1,1<<n)[None,:]>>np.arange(n)[:,None])&1)
        rows=[]
        for index,g in enumerate(atlas):
            if len(g)!=n-1:
                continue
            a=np.ones((n,n),dtype=np.int16)-np.eye(n,dtype=np.int16)
            for i,j in g.edges():
                a[i+1,j+1]=a[j+1,i+1]=-1
            energies=np.array([a[i,j] for i,j in pairs],dtype=np.int16)@products
            caps=[int((s*energies[full]).max()) for s in [1,-1]]
            certificates=[]
            for s in [1,-1]:
                values=np.full(1<<n,-n*n,dtype=np.int16)
                np.maximum.at(values,supports,s*energies)
                certificates.append(packing(values[1:],incidence))
            rows.append(dict(atlas_index=index, matrix=a.tolist(), P=caps[0],R=caps[1],
                             Q=max(caps),packing=certificates))
        minimum=min(row['Q'] for row in rows)
        minimizers=[r for r in rows if r['Q']==minimum]
        key=lambda r:max(F(c['value']) for c in r['packing'])
        overall=min(rows,key=key)
        selected=min(minimizers,key=key)
        record=dict(n=n,atlas_representatives=len(rows),M=minimum,
                    minimizing_representatives=len(minimizers),
                    minimum_packing=str(key(overall)),
                    minimum_packing_among_cap_minimizers=str(key(selected)),
                    best_packing=overall,best_cap_minimizing_packing=selected,
                    cap_minimizer_packing_values=sorted(set(str(key(r)) for r in minimizers)),
                    rows=rows,seconds=time.monotonic()-start)
        records.append(record)
        print(json.dumps({k:v for k,v in record.items() if k not in
                          ['rows','best_packing','best_cap_minimizing_packing']}),flush=True)
    return dict(status='integer caps and exact rational primal-dual certificates',
                coverage='NetworkX graph atlas, all first-row-positive signings up to residual vertex permutation',
                networkx_version=nx.__version__,atlas_counts=counts,records=records)


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    ap.add_argument('--max-n',type=int,default=8)
    ap.add_argument('--output',default='computations/results/flatify_director_small_packing_census_2026_09_07.json')
    args=ap.parse_args()
    assert 3<=args.max_n<=8
    data=run(args.max_n)
    data['source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    Path(args.output).write_text(json.dumps(data,indent=2)+'\n')
