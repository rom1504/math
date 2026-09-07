"""Independent integer reconstruction of saved actual-child bridge witnesses.

Includes the full parent histogram, child-energy shell envelope, and exact
one/two-edge bridge edit census. No numerical optimizer enters this check.
"""
import itertools
import json
from pathlib import Path
import numpy as np


def cube(n):
    return np.array([(1,)+v for v in itertools.product((-1,1),repeat=n-1)],dtype=np.int16)


def run():
    base=Path(__file__).resolve().parent/'results'
    source=json.loads((base/'flatify_adversary_2026_09_07_hadamard_bridge_phase_search.json').read_text())
    records=[]
    for row in source['records']:
        c=np.array(row['best']['parent_matrix'],dtype=np.int16); n=len(c)//2
        assert np.array_equal(c,c.T) and np.all(np.diag(c)==0)
        assert np.all(np.abs(c+np.eye(2*n,dtype=np.int16))==1)
        x=cube(n); y=cube(2*n); a=c[:n,:n]; d=c[n:,n:]; b=c[:n,n:]
        ha=np.sum((x@a)*x,axis=1)//2; hd=np.sum((x@d)*x,axis=1)//2
        internal=np.abs(ha[:,None]+hd[None,:]); bridge=x@b@x.T
        energy=np.sum((y@c)*y,axis=1)//2
        cap=int(np.max(np.abs(energy)))
        assert cap==int((internal+np.abs(bridge)).max())==row['best']['cap']
        assert np.array_equal(b@b.T,n*np.eye(n,dtype=np.int16))
        values,counts=np.unique(energy,return_counts=True)
        shells=[dict(child_sum=int(t),bridge_max=int(np.abs(bridge[internal==t]).max()),
                     pair_count=int(np.sum(internal==t)),
                     maximizing_pairs=int(np.sum((internal==t)&(internal+np.abs(bridge)==cap))))
                for t in np.unique(internal)]
        local={}
        if n==8:
            change=np.array([(-2*b[i,j]*np.outer(x[:,i],x[:,j])).ravel()
                             for i in range(n) for j in range(n)],dtype=np.int16)
            original=bridge.ravel(); fixed=internal.ravel()
            for k in (1,2):
                histogram={}
                for edges in itertools.combinations(range(n*n),k):
                    q=int(np.max(fixed+np.abs(original+np.sum(change[list(edges)],axis=0))))
                    histogram[q]=histogram.get(q,0)+1
                local[k]=histogram
        records.append(dict(child_order=n,child_caps=[int(np.abs(ha).max()),int(np.abs(hd).max())],
            parent_cap=cap,bridge_cap=int(np.abs(bridge).max()),parent=c.tolist(),bridge=b.tolist(),
            energy_histogram={int(v):int(q) for v,q in zip(values,counts)},
            shell_envelope=shells,exact_bridge_edit_cap_histograms=local))
    out=dict(status='PASS: INTEGER COMPLETE PROJECTIVE CUBES',records=records)
    (base/'flatify_adversary_2026_09_07_hadamard_bridge_witness_check.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in records[-1].items() if k not in ('parent','bridge','energy_histogram')}),flush=True)


if __name__=='__main__': run()
