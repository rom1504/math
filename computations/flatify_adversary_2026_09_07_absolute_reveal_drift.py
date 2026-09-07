"""Full-sign exhaustive numerical drift decomposition; not a proof certificate.

All signings modulo vertex switching are included. Successor profiles retain
the SAME vertex labels when evaluating a frozen current minimizing signing.
This is essential: independently canonicalizing successors would corrupt the
frozen-branch / optimizer-selection decomposition.
"""
import itertools
import json
import math
from pathlib import Path
import numpy as np


def signs(d):
    z=np.arange(1<<d,dtype=np.uint64)[:,None]
    return 1-2*((z>>np.arange(d,dtype=np.uint64))&1).astype(np.int8)


def setup(n):
    edges=list(itertools.combinations(range(n),2))
    x=np.column_stack((np.ones(1<<(n-1),dtype=np.int8),signs(n-1)))
    xx=np.array([x[:,i]*x[:,j] for i,j in edges],dtype=float)
    free=[e for e,(i,j) in enumerate(edges) if i]
    aa=np.ones((1<<len(free),len(edges)))
    aa[:,free]=signs(len(free))
    return edges,xx,aa


def profile(n,m,labels,edges):
    u=labels.count(0); r=labels.count(1); R=m-r
    a=(n-1)/(m-1); b=(n-1)/(n-m-1)
    out=[]
    for i,j in edges:
        s,t=labels[i],labels[j]
        if s and t:
            v=(a if s==1 else b) if s==t else 0
        elif s or t:
            v=a*R/u if s+t==1 else b*(u-R)/u
        else:
            v=(a*R*(R-1)+b*(u-R)*(u-R-1))/(u*(u-1))
        out.append(v)
    out=np.array(out)
    for i in range(n):
        assert abs(sum(out[e] for e,pair in enumerate(edges) if i in pair)-(n-1))<1e-10
    return out


def pressures(v,beta,n,xx,aa):
    energy=(aa*(beta*np.sqrt(v/(n-1))))@xx
    plus=np.logaddexp.reduce(energy,axis=1)+math.log(2)
    minus=np.logaddexp.reduce(-energy,axis=1)+math.log(2)
    absolute=np.logaddexp(plus,minus)-math.log(2)
    paired=(plus+minus)/2
    return absolute,paired


def run():
    results=[]
    for n in [4,6,7]:
        m=n//2
        edges,xx,aa=setup(n)
        for beta in [.7,1.5,3.]:
            records=[]
            for k in range(n):
                for r in range(max(0,k-(n-m)),min(k,m)+1):
                    labels=[1]*r+[-1]*(k-r)+[0]*(n-k)
                    v=profile(n,m,labels,edges)
                    fs,ws=pressures(v,beta,n,xx,aa)
                    index=int(fs.argmin()); f=float(fs[index])
                    chance=(m-r)/(n-k)
                    expected_min=0.; expected_frozen=0.
                    avg_v=np.zeros_like(v)
                    for label,p in [(1,chance),(-1,1-chance)]:
                        if p==0: continue
                        nxt=labels.copy(); nxt[k]=label
                        vv=profile(n,m,nxt,edges)
                        ff,_=pressures(vv,beta,n,xx,aa)
                        expected_min+=p*float(ff.min())
                        expected_frozen+=p*float(ff[index])
                        avg_v+=p*vv
                    assert np.max(np.abs(avg_v-v))<1e-10
                    selection=expected_frozen-expected_min
                    assert selection>-1e-10
                    probability=math.comb(m,r)*math.comb(n-m,k-r)/math.comb(n,k)
                    records.append(dict(k=k,r=r,probability=probability,
                        absolute_min=f,paired_min=float(ws.min()),
                        optimizing_signing_index=index,
                        optimized_drift=expected_min-f,
                        frozen_drift=expected_frozen-f,
                        selection_loss=selection))
            initial=records[0]['absolute_min']
            terminal,_=pressures(profile(n,m,[1]*m+[-1]*(n-m),edges),beta,n,xx,aa)
            endpoint=float(terminal.min())-initial
            telescoped=sum(z['probability']*z['optimized_drift'] for z in records)
            assert abs(endpoint-telescoped)<1e-9
            summary=dict(n=n,m=m,beta=beta,endpoint_drift=endpoint,
                total_frozen_drift=sum(z['probability']*z['frozen_drift'] for z in records),
                total_selection_loss=sum(z['probability']*z['selection_loss'] for z in records),
                min_conditional_drift=min(z['optimized_drift'] for z in records),
                max_conditional_drift=max(z['optimized_drift'] for z in records))
            print(json.dumps(summary),flush=True)
            results.append(dict(summary=summary,states=records))
    return dict(status='EXHAUSTIVE FULL SIGNINGS MODULO SWITCHING; FLOATING PRESSURES ONLY',records=results)


if __name__=='__main__':
    output=run()
    target=Path(__file__).resolve().parent/'results'/'flatify_adversary_2026_09_07_absolute_reveal_drift.json'
    target.write_text(json.dumps(output,indent=2)+'\n')
