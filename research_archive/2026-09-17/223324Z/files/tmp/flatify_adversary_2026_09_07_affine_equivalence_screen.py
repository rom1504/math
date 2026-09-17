"""Exact simultaneous rank-one tile test for affine-net realizability.

At dimension q^2, an affine representation is equivalent to row partitions
into q groups of q with every cross-group tile rank one. Enumerating a group
containing row zero is complete; projective column signatures force all other
partitions. Input bases are independently checked before the exact screen.
"""
import itertools
import json
from pathlib import Path
import sys
import time
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'computations'))
from flatify_construct_2026_09_07_kerdock_diagnostic import bases16
from flatify_construct_2026_09_07_mub_diagnostic import bases4

q=4;d=q*q

def check_bases(bases):
    for i,B in enumerate(bases):
        assert np.array_equal(B@B.T,d*np.eye(d,dtype=np.int64))
        for C in bases[:i]:assert np.all(np.abs(B@C.T)==q)

def partition_signatures(matrix):
    # Columns projectively equivalent; entries are literal +-1.
    signatures=matrix*matrix[:1,:]
    groups={}
    for j in range(d):groups.setdefault(tuple(signatures[:,j]),[]).append(j)
    groups=list(groups.values())
    return groups if len(groups)==q and all(len(g)==q for g in groups) else None

def rank1(tile):return np.array_equal(tile,tile[:,:1]*tile[:1,:]*tile[0,0])

def reconstruct_rotation(bases,parts):
    # U_i=bases_i/4; the common right rotation is W/8, all W entries integer.
    cross=(bases[0]@bases[1].T)//q
    columns=[]
    for P in parts[0]:
        for Q in parts[1]:
            signs=cross[P,Q[0]]
            columns.append(signs@bases[0][P,:])
    W=np.array(columns,dtype=np.int64).T
    assert np.array_equal(W.T@W,64*np.eye(d,dtype=np.int64))
    for B,partition in zip(bases,parts):
        transformed=B@W
        assert np.all(np.isin(transformed,[0,-16,16]))
        assert np.all(np.sum(transformed!=0,axis=1)==q)
        supports=[]
        for group in partition:
            assert all(np.array_equal(transformed[j]!=0,transformed[group[0]]!=0) for j in group)
            supports.append(set(np.flatnonzero(transformed[group[0]])))
        assert len(set.union(*supports))==d
    return W.tolist()

def find_partitions(bases):
    k=len(bases)
    crosses={(i,j):(bases[i]@bases[j].T)//q for i in range(k) for j in range(k) if i!=j}
    for rest in itertools.combinations(range(1,d),q-1):
        R=(0,)+rest
        parts=[None]
        for j in range(1,k):
            part=partition_signatures(crosses[0,j][list(R),:])
            if part is None:break
            parts.append(part)
        if len(parts)!=k:continue
        parts[0]=partition_signatures(crosses[0,1][:,parts[1][0]].T)
        if parts[0] is None:continue
        if all(rank1(crosses[i,j][np.ix_(P,Q)])
               for i in range(k) for j in range(i+1,k)
               for P in parts[i] for Q in parts[j]):
            return parts
    return None

start=time.monotonic()
affine=[2*V for V in bases4()]
check_bases(affine)
baseline=find_partitions(affine)
assert baseline is not None
baseline_rotation=reconstruct_rotation(affine,baseline)
kerdock=bases16();check_bases(kerdock)
cases=[]
for ids in itertools.combinations(range(len(kerdock)),5):
    witness=find_partitions([kerdock[i] for i in ids])
    rotation=reconstruct_rotation([kerdock[i] for i in ids],witness) if witness is not None else None
    cases.append(dict(ids=ids,affine_net_equivalent=witness is not None,partitions=witness,
                      common_rotation_scaled_by8=rotation))
six=[tuple(range(6))]
for ids in six:assert find_partitions([kerdock[i] for i in ids]) is None
out=dict(status='exact complete row-partition screen',dimension=d,q=q,
         baseline_affine_partitions=baseline,five_basis_cases=cases,
         affine_five_count=sum(r['affine_net_equivalent'] for r in cases),
         nonaffine_five_count=sum(not r['affine_net_equivalent'] for r in cases),
         six_basis_negative_controls=six,seconds=time.monotonic()-start)
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k not in ['baseline_affine_partitions','five_basis_cases']},indent=2))
