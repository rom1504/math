"""Exhaustive floating audit of the two-orientation cap-preserving reveal.

The internal revealed-positive edges get their own orientation variable.
All other edges share a second orientation. Actual signings are exhaustively
minimized modulo switching; successor labels are kept aligned for frozen
signing comparisons. No asymptotic sign claim follows from these tests.
"""
import json
import math
from pathlib import Path
import numpy as np
from flatify_adversary_2026_09_07_absolute_reveal_drift import setup,profile


def pressure(v,labels,beta,n,edges,xx,aa):
    mask=np.array([labels[i]==labels[j]==1 for i,j in edges])
    weighted=aa*(beta*np.sqrt(v/(n-1)))
    ep=(weighted*mask)@xx
    er=(weighted*(~mask))@xx
    terms=np.logaddexp(ep,-ep)+np.logaddexp(er,-er)-2*math.log(2)
    return np.logaddexp.reduce(terms,axis=1)+math.log(2)


def run():
    output=[]
    for n in [4,6,7]:
        m=n//2
        edges,xx,aa=setup(n)
        for beta in [.7,1.5,3.]:
            records=[]
            for k in range(n):
                for r in range(max(0,k-(n-m)),min(k,m)+1):
                    labels=[1]*r+[-1]*(k-r)+[0]*(n-k)
                    v=profile(n,m,labels,edges)
                    ff=pressure(v,labels,beta,n,edges,xx,aa)
                    index=int(ff.argmin()); f=float(ff[index])
                    chance=(m-r)/(n-k)
                    future=0.; frozen=0.
                    for label,p in [(1,chance),(-1,1-chance)]:
                        if not p: continue
                        nxt=labels.copy(); nxt[k]=label
                        vv=profile(n,m,nxt,edges)
                        gg=pressure(vv,nxt,beta,n,edges,xx,aa)
                        future+=p*float(gg.min()); frozen+=p*float(gg[index])
                    assert frozen>=future-1e-9
                    probability=math.comb(m,r)*math.comb(n-m,k-r)/math.comb(n,k)
                    records.append(dict(k=k,r=r,probability=probability,pressure=f,
                        optimized_drift=future-f,frozen_drift=frozen-f,
                        selection_loss=frozen-future,optimizing_signing_index=index))
            terminal_labels=[1]*m+[-1]*(n-m)
            terminal=float(pressure(profile(n,m,terminal_labels,edges),terminal_labels,beta,n,edges,xx,aa).min())
            child_total=0.
            for size in [m,n-m]:
                ce,cx,ca=setup(size)
                child_total+=float(pressure(np.ones(len(ce)),[0]*size,beta,size,ce,cx,ca).min())
            assert abs(terminal-child_total)<1e-9
            endpoint=terminal-records[0]['pressure']
            assert abs(endpoint-sum(z['probability']*z['optimized_drift'] for z in records))<1e-9
            summary=dict(n=n,m=m,beta=beta,endpoint_drift=endpoint,
                total_frozen_drift=sum(z['probability']*z['frozen_drift'] for z in records),
                total_selection_loss=sum(z['probability']*z['selection_loss'] for z in records),
                min_conditional_drift=min(z['optimized_drift'] for z in records),
                max_conditional_drift=max(z['optimized_drift'] for z in records),
                terminal_pressure=terminal,child_absolute_pressure_sum=child_total)
            print(json.dumps(summary),flush=True)
            output.append(dict(summary=summary,states=records))
    return dict(status='EXHAUSTIVE SIGNINGS; FLOATING PRESSURES ONLY',records=output)


if __name__=='__main__':
    output=run()
    target=Path(__file__).resolve().parent/'results'/'flatify_adversary_2026_09_07_grouped_reveal_drift.json'
    target.write_text(json.dumps(output,indent=2)+'\n')
