"""Bounded exact ground-field diagnostic on the four stored optimal A8 types."""
import itertools
import json
from pathlib import Path
import numpy as np

source=Path('computations/results/flatify_adversary_2026_09_07_hadamard_bridge_phase_search.json')
data=json.loads(source.read_text())
spins=np.array(list(itertools.product([-1,1],repeat=8)),dtype=np.int64)
records=[]
for item in data['records'][1]['child_representatives']:
    A=np.array(item['matrix'],dtype=np.int64)
    energies=np.einsum('bi,ij,bj->b',spins,A,spins)//2
    assert np.max(abs(energies))==10
    patterns=set()
    for x,e in zip(spins,energies):
        if abs(e)==10:
            patterns.add(tuple(sorted((int(np.sign(e))*x*(A@x)).tolist())))
    records.append({'atlas_index':item['atlas_index'],'patterns':sorted(patterns)})
result={'status':'EXACT FINITE DIAGNOSTIC; NO ASYMPTOTIC INFERENCE','records':records}
Path('computations/results/flatify_adversary_2026_09_07_small_ground_fields.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
