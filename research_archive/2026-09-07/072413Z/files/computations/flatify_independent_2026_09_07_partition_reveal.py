"""Exact-enumeration numerical pressure audit for a hidden partition reveal.

The variance profiles and state probabilities are exact rational quantities.
Pressure minimization exhausts switching-gauged signings; transcendental
pressure evaluations are floating point and are NOT proof certificates.
"""
from fractions import Fraction as F
import itertools
import json
import math
from pathlib import Path
import numpy as np


def profile(n,m,k,r):
    """r revealed positives, k-r negatives; target positive class has size m."""
    u=n-k
    a=F(n-1,m-1)
    b=F(n-1,n-m-1)
    remaining=m-r
    kinds=[1]*r+[-1]*(k-r)+[0]*u
    edges=list(itertools.combinations(range(n),2))
    v=[]
    for i,j in edges:
        x,y=kinds[i],kinds[j]
        if x and y:
            value=(a if x==1 else b) if x==y else F(0)
        elif x or y:
            known=x or y
            value=a*remaining/u if known==1 else b*(u-remaining)/u
        else:
            value=(a*remaining*(remaining-1)+b*(u-remaining)*(u-remaining-1))/(u*(u-1))
        v.append(value)
    for i in range(n):
        assert sum(v[e] for e,(j,l) in enumerate(edges) if i in (j,l))==n-1
    return edges,v


def signs(d):
    values=np.arange(2**d,dtype=np.uint64)[:,None]
    return (1-2*((values>>np.arange(d,dtype=np.uint64))&1).astype(np.int8))


def problem(n):
    edges=list(itertools.combinations(range(n),2))
    spins=np.column_stack([np.ones(2**(n-1),dtype=np.int8),signs(n-1)])
    spin_edges=np.asarray([spins[:,i]*spins[:,j] for i,j in edges],dtype=np.float64)
    free=[e for e,(i,j) in enumerate(edges) if i>0]
    patterns=np.ones((2**len(free),len(edges)),dtype=np.float64)
    patterns[:,free]=signs(len(free))
    return spin_edges,patterns


def pressure(v,beta,n,spin_edges,patterns):
    coupling=beta*np.sqrt(np.array(v,dtype=np.float64)/(n-1))
    energy=(patterns*coupling)@spin_edges
    answer=np.zeros(len(patterns))
    for orientation in [1,-1]:
        z=orientation*energy
        peak=z.max(axis=1)
        answer+=(peak+np.log(np.exp(z-peak[:,None]).sum(axis=1)))/2
    i=int(answer.argmin())
    return float(answer[i]+math.log(2)),i


def campaign():
    result=[]
    for n in [4,6,7]:
        m=n//2
        spin_edges,patterns=problem(n)
        for beta in [.2,.7,1.5,3.,6.]:
            states={}
            for k in range(n+1):
                for r in range(max(0,k-(n-m)),min(k,m)+1):
                    edges,v=profile(n,m,k,r)
                    f,best=pressure(v,beta,n,spin_edges,patterns)
                    states[k,r]={'pressure':f,'best_switching_code':best,
                                 'state_probability':str(F(math.comb(m,r)*math.comb(n-m,k-r),math.comb(n,k)))}
            rows=[]
            for k in range(n):
                average=0.
                average_drift=0.
                conditional=[]
                for r in range(max(0,k-(n-m)),min(k,m)+1):
                    source=states[k,r]
                    chance=F(m-r,n-k)
                    target=0.
                    if chance:
                        target+=float(chance)*states[k+1,r+1]['pressure']
                    if chance<1:
                        target+=float(1-chance)*states[k+1,r]['pressure']
                    drift=target-source['pressure']
                    weight=float(F(source['state_probability']))
                    average+=weight*source['pressure']
                    average_drift+=weight*drift
                    conditional.append({'r':r,'drift':drift,'chance_positive':str(chance)})
                rows.append({'k':k,'expected_pressure':average,'expected_drift':average_drift,
                             'conditional':conditional})
            total=sum(row['expected_drift'] for row in rows)
            endpoint=states[n,m]['pressure']-states[0,0]['pressure']
            assert abs(total-endpoint)<1e-10
            record={'n':n,'m':m,'beta':beta,'normalization':'sqrt(n-1)',
                    'endpoint_difference':endpoint,'rows':rows,
                    'states':{f'{k},{r}':v for (k,r),v in states.items()}}
            result.append(record)
            print(json.dumps({'n':n,'m':m,'beta':beta,'endpoint_difference':endpoint,
                              'drifts':[row['expected_drift'] for row in rows]}),flush=True)
    return {'status':'EXHAUSTIVE SIGNINGS, FLOATING PRESSURES; NOT ASYMPTOTIC',
            'records':result}


if __name__=='__main__':
    result=campaign()
    Path('computations/results/flatify_independent_2026_09_07_partition_reveal.json').write_text(
        json.dumps(result,indent=2)+'\n')
