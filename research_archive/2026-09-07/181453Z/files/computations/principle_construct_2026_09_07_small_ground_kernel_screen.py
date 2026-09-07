"""Bounded exact-seed screen; numerical eigenvalues are diagnostics only."""
import json
from pathlib import Path
import numpy as np

base=Path('computations/results')
paths=[(n,f'exact_m{n}.json','matrix') for n in range(3,11)]
paths += [(11,'nested_10_in_11_cap17.json','matrix'),
          (12,'extension_nested_m11_to_12.json','parent_matrix'),
          (13,'bridge_6_7_sign1_cap20.json','parent_matrix'),
          (14,'conference_completion_m13.json','conference_matrix')]
records=[]
for n,source,key in paths:
    a=np.array(json.loads((base/source).read_text())[key],dtype=np.int64)
    assert a.shape==(n,n) and np.array_equal(a,a.T) and not np.any(a.diagonal())
    spins=1-2*((np.arange(1<<(n-1))[:,None]>>np.arange(n))&1)
    energy=np.einsum('bi,ij,bj->b',spins,a,spins)//2
    for polarity in [1,-1]:
        aa=polarity*a; ee=polarity*energy; p=int(ee.max()); ground=spins[ee==p]
        gram=ground@aa@ground.T
        assert np.all(np.abs(gram)<=2*p)
        normalized=gram/(2*p)
        record={'n':n,'source':source,'polarity':polarity,'absolute_cap_exact':int(abs(energy).max()),
                'positive_cap_exact':p,'projective_ground_count':len(ground),
                'normalized_gram_min_eigenvalue_DISPLAY':float(np.linalg.eigvalsh(normalized)[0]),
                'cosh_kernel_min_eigenvalues_DISPLAY':{str(t):float(np.linalg.eigvalsh(np.cosh(t*normalized))[0]) for t in [.1,1.,4.]}}
        witness=None
        for center in range(len(ground)):
            near=np.flatnonzero(np.abs(gram[center])==2*p)
            for i in near:
                other=near[np.abs(gram[i,near])<2*p]
                if len(other):
                    j=int(other[0]); ids=[int(i),j,center]
                    witness={'indices':ids,'spins':ground[ids].tolist(),'integer_gram':gram[np.ix_(ids,ids)].tolist(),
                             'weights':[1,1,-2],'certificate':'2*(cosh(t*abs(G12)/(2P))-cosh(t))<0 for every t>0'}
                    break
            if witness: break
        record['exact_three_point_saturated_witness']=witness
        records.append(record)
print(json.dumps({'status':'FINITE EXACT INPUTS; EIGENVALUES NUMERICAL; SATURATED TRIPLES EXACT; NO ASYMPTOTIC INFERENCE',
                  'records':records},indent=2))
