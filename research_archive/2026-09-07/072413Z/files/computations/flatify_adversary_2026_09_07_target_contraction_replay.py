"""Independent numerical stress replay of target-contraction sign recovery.

Checks signed/high-leverage/full-rank perturbations, mean coefficient
feasibility, full-cube contraction, variance/edit bounds, and actual signs.
Floating validation, not a replacement for the exact algebraic proof.
"""
import json
import math
from pathlib import Path
import numpy as np


def run():
    rng=np.random.default_rng(9070704)
    records=[]
    for n in [4,7,10,12]:
        ei,ej=np.triu_indices(n,1)
        z=np.arange(1<<(n-1),dtype=np.uint64)[:,None]
        spins=np.column_stack((np.ones(len(z)),1-2*((z>>np.arange(n-1,dtype=np.uint64))&1).astype(np.int8)))
        products=(spins[:,ei]*spins[:,ej]).T
        for case in range(12):
            av=rng.choice([-1,1],len(ei))
            A=np.zeros((n,n)); A[ei,ej]=av; A[ej,ei]=av
            if case%3==0:
                u=rng.normal(size=n); u[0]*=20
                u/=np.linalg.norm(u)
                delta=(-1 if case%2 else 1)*(case+1)*np.outer(u,u)
            else:
                U,_=np.linalg.qr(rng.normal(size=(n,n)))
                ev=rng.uniform(-1,1,n)*(case+1)/3
                if case%3==1: ev[2:]=0
                delta=(U*ev)@U.T
            eig,U=np.linalg.eigh(delta)
            absdelta=(U*np.abs(eig))@U.T
            ell=np.maximum(np.diag(absdelta),0)
            L=float(ell.sum())
            s=1/np.sqrt(1+ell)
            target=av-delta[ei,ej]
            mean=s[ei]*s[ej]*target
            assert np.max(np.abs(mean))<=1+1e-10
            cap_target=float(np.max(np.abs(target@products)))
            cap_mean=float(np.max(np.abs(mean@products)))
            assert cap_mean<=cap_target+1e-9
            variance=float(np.sum(1-mean**2))
            edits=float(np.sum((1-av*mean)/2))
            assert variance<=2*(n-1)*L+1e-8
            assert edits<=(n-1)*L/2+1e-8
            a=(n+2)*math.log(2)
            error=2*math.sqrt(n*L*a)+4*a/3
            witness=None
            for attempt in range(100):
                realized=np.where(rng.random(len(ei))<(1+mean)/2,1,-1)
                discrepancy=float(np.max(np.abs((realized-mean)@products)))
                changed=int(np.sum(realized!=av))
                if discrepancy<=error and changed<=2*n*L:
                    witness=realized
                    break
            assert witness is not None
            records.append(dict(n=n,case=case,nuclear_budget=L,
                input_signs=av.tolist(),perturbation=delta.tolist(),mean_edges=mean.tolist(),
                output_signs=witness.tolist(),cap_target=cap_target,cap_mean=cap_mean,
                cap_output=float(np.max(np.abs(witness@products))),
                variance=variance,variance_upper=2*(n-1)*L,
                expected_edits=edits,actual_edits=changed,
                centered_cap=discrepancy,rounding_bound=error,attempts=attempt+1))
    summary=dict(status='PASS: FLOATING ALGEBRA/STRESS REPLAY, NOT A PROOF CERTIFICATE',
        cases=len(records),largest_mean_amplitude=max(max(abs(t) for t in x['mean_edges']) for x in records),
        largest_contraction_ratio=max(x['cap_mean']/x['cap_target'] for x in records))
    print(json.dumps(summary),flush=True)
    return dict(summary=summary,records=records)


if __name__=='__main__':
    result=run()
    target=Path(__file__).resolve().parent/'results'/'flatify_adversary_2026_09_07_target_contraction_replay.json'
    target.write_text(json.dumps(result,indent=2)+'\n')
