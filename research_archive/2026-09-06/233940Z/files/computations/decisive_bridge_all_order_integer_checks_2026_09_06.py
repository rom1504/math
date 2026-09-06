"""Exact integer replay of H12 and small restricted-weave identities."""
import json
from pathlib import Path
import numpy as np
from continued_convergence_restricted_weave_2026_09_06 import instance


def main():
    q=11
    squares={i*i%q for i in range(1,q)}
    character=np.array([0]+[1 if i in squares else -1 for i in range(1,q)],dtype=np.int64)
    c=character[(np.arange(q)[:,None]-np.arange(q)[None,:])%q]
    r=np.zeros((q+1,q+1),dtype=np.int64)
    r[0,1:]=1;r[1:,0]=-1;r[1:,1:]=c
    h12=np.eye(q+1,dtype=np.int64)+r
    assert np.all(abs(h12)==1)
    assert np.array_equal(h12@h12.T,12*np.eye(12,dtype=np.int64))
    a,b,t,s,perms,rows,_=instance(4,3,20260906)
    count=0
    for bits in range(1<<12):
        x=np.array([1-2*((bits>>j)&1) for j in range(12)],dtype=np.int64)
        h=np.einsum("iaj,ia->ij",t,x.reshape(4,3))
        energy=int(x@b@x)
        assert energy==int(np.sum(s*h*h.T))
        assert int(np.sum(h*h))==4*4*3
        for sigma in (-1,1):
            assert int(np.sum((h-sigma*s*h.T)**2))==2*(4*4*3-sigma*energy)
        assert abs(int(x@a@x)-energy)<=12
        count+=1
    result={"status":"exact integer identities, not an asymptotic cap certificate",
            "H12":h12.tolist(),"paley_order":12,"m":4,"retained_per_fibre":3,
            "seed":20260906,"spin_assignments_checked":count,
            "row_selectors":rows.tolist(),"column_permutations":perms.tolist(),"outer_signs":s.tolist()}
    Path("computations/decisive_bridge_all_order_integer_checks_2026_09_06.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k not in ("H12","row_selectors","column_permutations","outer_signs")},indent=2))


if __name__=="__main__":main()
