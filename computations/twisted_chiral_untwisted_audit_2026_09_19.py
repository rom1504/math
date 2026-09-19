#!/usr/bin/env python3
"""Independent direct-spin comparison B=A versus all certified twists."""
import json
from pathlib import Path
import numpy as np

import twisted_chiral_driver_2026_09_18 as driver
from twisted_chiral_summarize_2026_09_18 import collect,FAMILY_TAGS


def audit(a):
    a=np.asarray(a,dtype=np.int64);r=len(a)
    ym=np.arange(1<<r,dtype=np.int64)
    xm=np.arange(1<<(r-1),dtype=np.int64)<<1
    y=1-2*((ym[:,None]>>np.arange(r))&1)
    x=1-2*((xm[:,None]>>np.arange(r))&1)
    qx=np.einsum('bi,ij,bj->b',x,a,x)//2
    qy=np.einsum('bi,ij,bj->b',y,a,y)//2
    base=qx[:,None]-qy[None,:]+x@a@y.T
    masks=xm[:,None]^ym[None,:]
    lo=np.full(1<<r,10000,dtype=np.int64)
    hi=np.full(1<<r,-10000,dtype=np.int64)
    np.minimum.at(lo,masks.ravel(),base.ravel())
    np.maximum.at(hi,masks.ravel(),base.ravel())
    # Rows range over every d, columns over every sector t; no chiral quotient.
    shifts=y@y.T
    caps=np.maximum(np.max(np.abs(lo[None,:]+shifts),axis=1),
                    np.max(np.abs(hi[None,:]+shifts),axis=1))
    best=int(np.argmin(caps))
    return dict(cap=int(caps[best]),p=list(range(r)),s=[1]*r,d=y[best].tolist(),
        exact_matching_count=1<<r,projective_parent_state_count=1<<(2*r-1),
        minimizing_matching_count=int(np.count_nonzero(caps==caps[best])),
        untwisted_core_cap=int(np.max(np.abs(base))),
        untwisted_conditional_radius=int(np.max(hi-lo)//2),
        all_plus_matching_cap=int(caps[0]))


def main():
    driver.ensure_binaries();records=[]
    families=collect(FAMILY_TAGS,'complete_family')
    extra=json.loads((driver.ROOT/'computations/results/twisted_chiral_2026_09_18_alternative_conference.json').read_text())['records'][0]
    for label,f in [(f'order{r}_class{c}',v) for (r,c),v in sorted(families.items())]+[('order10_nonoptimal_conference',extra)]:
        row=audit(f['child_matrix']);driver.verify(f['child_matrix'],row)
        row.update(label=label,full_twisted_family_minimum=f['cap'],
            twist_improvement=row['cap']-f['cap'],
            evidence='independent direct integer matrix-product audit over x0=+1 and unrestricted y; every d tested')
        records.append(row)
        print(label,row['child_cap'],row['cap'],f['cap'],row['all_plus_matching_cap'],flush=True)
    output=driver.ROOT/'computations/results/twisted_chiral_untwisted_audit_2026_09_19.json'
    output.write_text(json.dumps({'records':records},indent=2)+'\n')


if __name__=='__main__':main()
