"""Exhaustive small-order regression for the balanced-scar identities."""
import itertools
import json
from pathlib import Path
import numpy as np

root = Path(__file__).resolve().parents[1]
m, q, ell = 4, 4, 1
n = m*q
rng = np.random.default_rng(202609072030)
balanced = np.asarray([v for v in itertools.product((-1, 1), repeat=q)
                       if sum(v) == 0], dtype=np.int64)
words = np.asarray([np.concatenate(v) for v in
                    itertools.product(balanced, repeat=m)], dtype=np.int64)
base = np.zeros((n,n), dtype=np.int64)
factors = {}
for i in range(m):
    for j in range(i+1,m):
        left = balanced[rng.integers(len(balanced))].copy()
        right = balanced[rng.integers(len(balanced))].copy()
        factors[i,j] = left, right
        block = np.outer(left, right)
        base[i*q:(i+1)*q,j*q:(j+1)*q] = block
        base[j*q:(j+1)*q,i*q:(i+1)*q] = block.T
base_energies = np.einsum('bi,ij,bj->b', words, base, words)//2
base_cap = int(np.max(np.abs(base_energies)))
ferro = np.kron(np.eye(m,dtype=np.int64),
                np.ones((q,q),dtype=np.int64)-np.eye(q,dtype=np.int64))
ferro_ground = m*q*(q-1)//2
records = []
for selected in itertools.product(range(q),repeat=m):
    scar = base.copy()
    switches = 0
    target_change = 0
    for (i,j),(left,right) in factors.items():
        a,c = selected[i], selected[j]
        original = int(base[i*q+a,j*q+c])
        target_change += 1-original
        if original == 1:
            continue
        u = next(v for v in range(q) if v != a and left[v] == -left[a])
        v = next(w for w in range(q) if w != c and right[w] == -right[c])
        rows, cols = [i*q+a,i*q+u], [j*q+c,j*q+v]
        assert np.array_equal(scar[np.ix_(rows,cols)], [[-1,1],[1,-1]])
        for ii in rows:
            for jj in cols:
                scar[ii,jj] *= -1
                scar[jj,ii] = scar[ii,jj]
        switches += 1
    for i in range(m):
        for j in range(m):
            block = scar[i*q:(i+1)*q,j*q:(j+1)*q]
            assert np.all(block.sum(axis=0)==0)
            assert np.all(block.sum(axis=1)==0)
            if i != j:
                assert np.all(np.abs(block)==1)
    assert np.array_equal(scar,scar.T)
    z = np.zeros(n,dtype=np.int64)
    z[np.arange(m)*q+np.asarray(selected)] = 1
    x = 1-2*z
    old_z = int(z@base@z//2)
    new_z = int(z@scar@z//2)
    new_x = int(x@scar@x//2)
    assert new_z == old_z+target_change
    assert new_x == 4*new_z
    variance = sum(int(np.sum((x[i*q:(i+1)*q]-x[i*q:(i+1)*q].mean())**2))
                   for i in range(m))
    assert variance == 4*m*ell*(1-ell/q)
    scar_energies = np.einsum('bi,ij,bj->b', words, scar, words)//2
    cap = int(np.max(np.abs(scar_energies)))
    max_word_change = int(np.max(np.abs(scar_energies-base_energies)))
    assert max_word_change <= 8*switches
    assert abs(cap-base_cap) <= 8*ell*m*(m-1)//2
    assert np.count_nonzero(np.triu(scar-base,1)) == 4*switches
    parent_energy = int(x@(scar+ferro)@x//2)
    gain = parent_energy-ferro_ground
    assert gain == 4*new_z-2*m*ell*(q-ell)
    assert int(np.ones(n,dtype=np.int64)@(scar+ferro)@np.ones(n,dtype=np.int64)//2) == ferro_ground
    records.append(dict(selected=list(selected), switches=switches,
                        original_selected_energy=old_z, scar_selected_energy=new_z,
                        target_change=target_change, biased_bulk_energy=new_x,
                        centered_variance=variance, balanced_cap=cap,
                        max_balanced_word_energy_change=max_word_change,
                        ferro_completion_gain=gain))

assert sum(r['original_selected_energy'] for r in records) == 0
assert sum(r['scar_selected_energy'] for r in records) == len(records)*ell*m*(m-1)//2
assert sum(r['ferro_completion_gain'] for r in records) == len(records)*(4*ell*m*(m-1)//2-2*m*ell*(q-ell))
result = dict(status='PASS', m=m,q=q,ell=ell,
              selected_choices=len(records), balanced_words=len(words),
              base_matrix=base.tolist(), base_balanced_cap=base_cap,
              mean_scar_selected_energy=sum(r['scar_selected_energy'] for r in records)/len(records),
              balanced_cap_range=[min(r['balanced_cap'] for r in records),max(r['balanced_cap'] for r in records)],
              ferro_gain_range=[min(r['ferro_completion_gain'] for r in records),max(r['ferro_completion_gain'] for r in records)],
              records=records,
              scope='Exact finite regression only; the asymptotic scar theorem is proved analytically.')
target = root/'computations/results/principle_construct_2026_09_07_balanced_scar_check.json'
target.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in ('records','base_matrix')}))
