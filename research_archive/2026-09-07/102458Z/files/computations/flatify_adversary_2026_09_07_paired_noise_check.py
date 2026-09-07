"""Exact actual-child paired-noise sanity check and directed entropy display."""
from pathlib import Path
from fractions import Fraction as Q
import itertools
import json
import numpy as np
from scipy.linalg import hadamard
import mpmath as mp
from flatify_independent_2026_09_07_ternary_interval import endpoints, rational_iv

data=json.loads(Path('computations/results/flatify_adversary_2026_09_07_hadamard_bridge_phase_search.json').read_text())
A=np.array(data['records'][1]['child_representatives'][0]['matrix'],dtype=np.int64)
spins=np.array(list(itertools.product([-1,1],repeat=8)),dtype=np.int64)
energies=np.einsum('bi,ij,bj->b',spins,A,spins)//2
assert np.max(abs(energies))==10
index=int(np.argmax(abs(energies)))
x0=spins[index]
if energies[index]<0: A=-A
B=x0[:,None]*A*x0[None,:]
assert int(np.sum(B)//2)==10 and np.all(B.sum(axis=1)>=0)
groups=[[0,1],[2,3],[4],[5],[6],[7]]
delta=Q(1,10); alpha=1-2*delta
mean=Q(0); second=Q(0); count=0
H=hadamard(4).astype(np.int64)
for signs in itertools.product([-1,1],repeat=len(groups)):
    eta=np.ones(8,dtype=np.int64)
    probability=Q(1)
    for group,sign in zip(groups,signs):
        eta[group]=sign
        probability*=delta if sign==-1 else 1-delta
    spectrum=H.T@eta[:4]
    assert np.all(spectrum[[1,3]]==0)
    value=int(eta@B@eta//2)
    mean+=probability*value
    second+=probability*value*value
    count+=1
prediction=alpha**2*10+(1-alpha**2)*int(B[0,1]+B[2,3])
assert mean==prediction
mp.iv.dps=45
d=rational_iv(delta)
h=-d*mp.iv.ln(d)-(1-d)*mp.iv.ln(1-d)
entropy=rational_iv(Q(39,20))*h
assert endpoints(entropy)[0]>Q(633,1000)
result={'status':'PASS: EXACT ACTUAL CHILD AND DIRECTED ENTROPY',
        'child':A.tolist(),'ground_center':x0.tolist(),'child_cap':10,
        'noise_patterns':count,'exact_mean':str(mean),
        'predicted_mean':str(prediction),'exact_variance':str(second-mean*mean),
        'pair_entropy_interval':list(map(str,endpoints(entropy))),
        'asymptotic_defect_squared_lower':'9/32000'}
Path('computations/results/flatify_adversary_2026_09_07_paired_noise_check.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ['child','ground_center']},indent=2))
