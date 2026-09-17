"""Exact character-variance replay for independent rank-one dictionaries."""
import itertools
import json
import numpy as np

rng=np.random.default_rng(20260917)
reports=[]
for dimensions in ([(2,2)],[(2,3)],[(3,3)],[(2,2),(1,2)],[(2,3),(2,2)]):
    bit_count=sum(k+p for k,p in dimensions)
    underlying=np.asarray(list(itertools.product((-1,1),repeat=bit_count)),dtype=np.int64)
    chunks=[]
    characters=[]
    offset=0
    for k,p in dimensions:
        chunks.append((underlying[:,offset:offset+k,None]*underlying[:,None,offset+k:offset+k+p]).reshape(len(underlying),k*p))
        characters.extend((1<<(offset+i))^(1<<(offset+k+j)) for i in range(k) for j in range(p))
        offset+=k+p
    words=np.concatenate(chunks,axis=1)
    n=words.shape[1]
    edges=list(itertools.combinations(range(n),2))
    edge_features=np.asarray([words[:,i]*words[:,j] for i,j in edges]).T
    sizes=[k*p for k,p in dimensions]
    numerator=n*n+3*sum(s*s for s in sizes)-3*sum(k*p*(k+p) for k,p in dimensions)+2*n
    assert numerator%2==0
    sharp=numerator//2
    for trial in range(31):
        coefficients=np.ones(len(edges),dtype=np.int64) if trial==0 else rng.choice((-1,1),size=len(edges))
        energy=edge_features@coefficients
        class_sums={}
        for (i,j),coefficient in zip(edges,coefficients):
            character=characters[i]^characters[j]
            assert character
            class_sums[character]=class_sums.get(character,0)+int(coefficient)
        exact_variance=sum(value*value for value in class_sums.values())
        assert int(energy.sum())==0
        assert int(energy@energy)==len(words)*exact_variance
        assert exact_variance<=sharp
        if trial==0:
            assert exact_variance==sharp
    reports.append({"dimensions":dimensions,"order":n,"underlying_words":len(words),
                    "matrices_checked":31,"sharp_maximum_variance":sharp})
print(json.dumps({"status":"PASS exact rank-one energy variance","cases":reports},indent=2))
