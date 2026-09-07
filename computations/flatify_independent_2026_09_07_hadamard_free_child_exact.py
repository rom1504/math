"""Solver-free exhaustive replay of the order-8 Hadamard amalgamation test."""
import itertools
import json
import time
from pathlib import Path

import numpy as np


def direct_frames(X):
    gram = X @ X.T
    neighbors = [sum(1 << int(j) for j in np.flatnonzero(row == 0)) for row in gram]
    answer = []

    def visit(chosen, candidates):
        need = 8-len(chosen)
        if need == 0:
            answer.append(tuple(chosen))
            return
        while bin(candidates).count('1') >= need:
            bit = candidates & -candidates
            i = bit.bit_length()-1
            candidates ^= bit
            visit(chosen+[i],candidates & neighbors[i])

    visit([], (1 << 128)-1)
    return answer


def main():
    begin=time.time()
    X=np.array([[1]+[1 if z>>i&1 else -1 for i in range(7)] for z in range(128)],dtype=np.int64)
    edges=list(itertools.combinations(range(8),2))
    cuts=np.array([sum(1<<e for e,(i,j) in enumerate(edges) if x[i]!=x[j]) for x in X],dtype=np.uint32)
    table=np.array([bin(i).count('1') for i in range(65536)],dtype=np.uint8)

    def popcount(a):
        return table[a & 65535]+table[a >> 16]

    normalized=np.arange(1<<21,dtype=np.uint32) << 7
    for cut in cuts:
        count=popcount(normalized ^ cut)
        normalized=normalized[(count>=9)&(count<=19)]
    assert len(normalized)>0
    cap8=normalized.copy()
    for cut in cuts:
        count=popcount(cap8 ^ cut)
        cap8=cap8[(count>=10)&(count<=18)]
    assert len(cap8)==0
    print(json.dumps({'normalized_optimal_children':len(normalized),'cap8_children':len(cap8)}),flush=True)

    data=json.loads(Path('computations/results/flatify_adversary_2026_09_07_hadamard_bridge_phase_search.json').read_text())
    reps=data['records'][1]['child_representatives']
    orbit=set()
    for rep in reps:
        A=np.array(rep['matrix'],dtype=np.int64)
        assert np.all(A[0,1:]==1)
        for p in itertools.permutations(range(1,8)):
            p=(0,)+p
            orbit.add(sum(1<<e for e,(i,j) in enumerate(edges) if A[p[i],p[j]]<0))
    assert orbit == set(map(int,normalized)), 'atlas representatives do not exhaust the directly checked normalized minima'
    full=np.concatenate([normalized ^ c for c in cuts])
    assert len(np.unique(full))==len(full)
    allframes=direct_frames(X)
    assert len(allframes)==480 and len(set(allframes))==480
    records=[]
    conjugation_obstructions=[]
    for rep in reps:
        A=np.array(rep['matrix'],dtype=np.int64)
        hA=np.einsum('bi,ij,bj->b',X,A,X)//2
        energies=[sum(int(hA[i])**2 for i in frame) for frame in allframes]
        zero_counts=[sum(hA[i]==0 for i in frame) for frame in allframes]
        assert min(energies)>0
        conjugation_obstructions.append({'atlas_index':rep['atlas_index'],'projective_zero_energy_spins':int(np.sum(hA==0)),'maximum_zero_energy_columns_in_hadamard_frame':int(max(zero_counts)),'minimum_sum_column_energy_squares':int(min(energies))})
        for fi,frame in enumerate(allframes):
            B=X[list(frame)].T
            assert np.array_equal(B @ B.T,8*np.eye(8,dtype=int))
            bridge=np.abs(X @ B @ X.T)
            lower=np.maximum(np.max(-28+bridge-hA[:,None],axis=0),-10)
            upper=np.minimum(np.min(28-bridge-hA[:,None],axis=0),10)
            candidates=full
            checked=0
            for yi in np.argsort(upper-lower,kind='stable'):
                counts=popcount(candidates ^ cuts[yi]).astype(np.int16)
                h=28-2*counts
                candidates=candidates[(h>=lower[yi]) & (h<=upper[yi])]
                checked+=1
                if len(candidates)==0:
                    break
            assert len(candidates)==0, (rep['atlas_index'],fi,candidates.tolist())
            records.append({'atlas_index':rep['atlas_index'],'frame':list(frame),'shell_rows_checked':checked,'surviving_optimal_children':0})
        print(json.dumps({'completed_atlas':rep['atlas_index'],'frames':len(allframes),'elapsed':time.time()-begin}),flush=True)
    result={'status':'PASS: SOLVER-FREE EXHAUSTIVE INTEGER CHECK','normalized_optimal_children':len(normalized),'all_optimal_children':len(full),'minimum_child_cap':10,'projective_hadamard_frames':len(allframes),'tested_representative_frame_pairs':len(records),'cap28_completions':0,'conjugation_obstructions':conjugation_obstructions,'elapsed_seconds':time.time()-begin,'records':records}
    output=Path('computations/results/flatify_independent_2026_09_07_hadamard_free_child_exact.json')
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='records'}),flush=True)


if __name__=='__main__':
    main()
