"""Exact small-rank replay of the projective-character Gram identity."""
import itertools
import json
import numpy as np

records=[]
for rank in (1,2,3,4):
    reps=[]
    for x in itertools.product(range(7),repeat=rank):
        nonzero=next((a for a in x if a),None)
        if nonzero==1:
            reps.append(x)
    a=np.asarray(reps,dtype=np.int64)
    chi=np.array([0,1,1,-1,1,-1,-1],dtype=np.int64)
    c=chi[(a@a.T)%7]
    scale=7**(rank-1)
    assert len(reps)==(7**rank-1)//6
    assert np.array_equal(c@c.T,scale*np.eye(len(reps),dtype=np.int64))
    # Exact inversion for a nonconstant integer target phase vector.
    target=np.arange(len(reps),dtype=np.int64)**2-3
    numerator=c.T@target
    assert np.array_equal(c@numerator,scale*target)
    records.append(dict(rank=rank,projective_lines=len(reps),gram_scale=scale,passed=True))
print(json.dumps(records,indent=2))
