"""Exact finite dense bridge search with actual minimizing children.

Exhausts all row/column sign phases for each tested Hadamard bridge.
The permutation sample is bounded, not exhaustive. Child cap values and
all resulting parent caps are integer full-cube calculations.
"""
import json
import math
from pathlib import Path
import numpy as np
import networkx as nx


def spins(n):
    z=np.arange(1<<(n-1),dtype=np.uint64)[:,None]
    return np.column_stack((np.ones(len(z),dtype=np.int16),1-2*((z>>np.arange(n-1,dtype=np.uint64))&1).astype(np.int16)))


def hadamard(n):
    h=np.ones((1,1),dtype=np.int16)
    while len(h)<n: h=np.block([[h,h],[h,-h]])
    assert np.array_equal(h@h.T,n*np.eye(n,dtype=np.int16))
    return h


def children(n,target,x):
    result=[]
    for index,g in enumerate(nx.graph_atlas_g()):
        if len(g)!=n-1: continue
        a=np.ones((n,n),dtype=np.int16)-np.eye(n,dtype=np.int16)
        for i,j in g.edges: a[i+1,j+1]=a[j+1,i+1]=-1
        energies=np.sum((x@a)*x,axis=1)//2
        if int(np.max(np.abs(energies)))==target:
            result.append((index,a,energies))
    return result


def phase_search(h,ea,ed,x):
    k=len(x)
    bridge=np.abs(x@h@x.T)
    ix=np.arange(k)
    permutations=ix[:,None]^ix[None,:]
    eds=ed[permutations]
    best=None
    histogram={}
    for polarity in [1,-1]:
        for p in range(k):
            scores=np.abs(ea[ix^p][None,:,None]+polarity*eds[:,None,:])+bridge[None,:,:]
            caps=np.max(scores,axis=(1,2))
            values,counts=np.unique(caps,return_counts=True)
            for v,c in zip(values,counts): histogram[int(v)]=histogram.get(int(v),0)+int(c)
            q=int(caps.argmin()); cap=int(caps[q])
            if best is None or cap<best['cap']:
                best=dict(cap=cap,row_phase=p,column_phase=q,child_polarity=polarity)
    best['phase_cap_histogram']=histogram
    return best


def run():
    rng=np.random.default_rng(9070737)
    output=[]
    for n,target in [(4,4),(8,10)]:
        x=spins(n); h=hadamard(n); cs=children(n,target,x)
        records=[]; best=None
        tasks=[(i,j,np.arange(n),np.arange(n),'all child pairs, identity permutation')
               for i in range(len(cs)) for j in range(len(cs))]
        for attempt in range(12):
            tasks.append((int(rng.integers(len(cs))),int(rng.integers(len(cs))),
                          rng.permutation(n),rng.permutation(n),'sampled permutations'))
        for i,j,rp,cp,kind in tasks:
            hp=h[rp][:,cp]
            result=phase_search(hp,cs[i][2],cs[j][2],x)
            result.update(child_indices=[i,j],atlas_indices=[cs[i][0],cs[j][0]],
                          row_permutation=rp.tolist(),column_permutation=cp.tolist(),kind=kind)
            records.append(result)
            if best is None or result['cap']<best['cap']:
                best=result.copy()
                b=x[result['row_phase']][:,None]*hp*x[result['column_phase']][None,:]
                parent=np.block([[cs[i][1],b],[b.T,result['child_polarity']*cs[j][1]]])
                xp=spins(2*n)
                direct=int(np.max(np.abs(np.sum((xp@parent)*xp,axis=1)//2)))
                assert direct==best['cap']
                best['parent_matrix']=parent.tolist()
                best['bridge_matrix']=b.tolist()
            print(json.dumps(dict(child_order=n,child_pair=[i,j],kind=kind,
                                  cap=result['cap'],best_so_far=best['cap'])),flush=True)
        output.append(dict(child_order=n,actual_child_cap=target,
            child_representatives=[dict(atlas_index=i,matrix=a.tolist()) for i,a,e in cs],
            favorable_weighted_target=2*math.sqrt((2*n-1)/(n-1))*target,
            best=best,records=records))
    return dict(status='EXACT INTEGER PHASE SEARCH; BOUNDED PERMUTATION SAMPLE; NO ASYMPTOTIC CLAIM',records=output)


if __name__=='__main__':
    result=run()
    target=Path(__file__).resolve().parent/'results'/'flatify_adversary_2026_09_07_hadamard_bridge_phase_search.json'
    target.write_text(json.dumps(result,indent=2)+'\n')
