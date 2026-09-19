#!/usr/bin/env python3
"""Exact switching-only minima and independent partial-reversal identity."""
import hashlib
import json
import subprocess
import numpy as np

import twisted_chiral_driver_2026_09_18 as base
from twisted_chiral_summarize_2026_09_18 import collect,FAMILY_TAGS


def lift(a,s,d):
    b=a*s[:,None]*s[None,:]+np.diag(d)
    return np.block([[a,b],[b,-a]])


def rotation_identity(a,s,d):
    n=len(a);j=s==-1
    partial=a.copy();partial[np.ix_(j,j)]*=-1
    rotation=np.eye(2*n,dtype=np.int64)
    for i in range(n):
        if j[i]:
            rotation[i,i]=rotation[i+n,i+n]=0
            rotation[i,i+n]=1;rotation[i+n,i]=-1
    assert np.array_equal(rotation.T@lift(a,s,d)@rotation,
                          lift(partial,np.ones(n,dtype=np.int64),d*s))
    return partial


def partial_reversal_caps(a):
    n=len(a);masks=np.arange(1<<n,dtype=np.int64)
    bits=(masks[:,None]>>np.arange(n))&1
    x=1-2*((np.arange(1<<(n-1),dtype=np.int64)[:,None]>>np.arange(n))&1)
    # This alternative projective chart fixes the LAST spin positive.
    edges=[(i,j) for i in range(n) for j in range(i+1,n)]
    products=np.asarray([x[:,i]*x[:,j] for i,j in edges])
    reversed_edges=np.asarray([a[i,j]*(1-2*(bits[:,i]*bits[:,j])) for i,j in edges]).T
    energies=reversed_edges@products
    caps=np.max(np.abs(energies),axis=1)
    histogram={str(int(c)):int(np.count_nonzero(caps==c)) for c in sorted(set(caps))}
    assert np.array_equal(caps,caps[::-1])
    return dict(reversal_count=1<<n,minimum_child_cap=int(caps.min()),
        minimizing_reversal_masks=np.flatnonzero(caps==caps.min()).tolist(),cap_histogram=histogram)


def main():
    base.ensure_binaries()
    source=base.ROOT/'computations/twisted_chiral_switching_only_2026_09_19.cpp'
    dependency=base.ROOT/'computations/twisted_chiral_search_frozen_2026_09_18.cpp'
    binary=base.SCRATCH/'switching_only_2026_09_19'
    subprocess.run(['g++','-O3','-march=native','-std=c++17',str(source),'-o',str(binary)],check=True)
    families=collect(FAMILY_TAGS,'complete_family')
    conference=json.loads((base.ROOT/'computations/results/twisted_chiral_2026_09_18_alternative_conference.json').read_text())['records'][0]
    records=[]
    for label,f in [(f'order{r}_class{c}',v) for (r,c),v in sorted(families.items())]+[('order10_nonoptimal_conference',conference)]:
        a=np.asarray(f['child_matrix'],dtype=np.int64)
        output=subprocess.check_output([str(binary)],input=base.matrix_input(a),text=True)
        items=[json.loads(line) for line in output.splitlines()]
        record=dict(next(row for row in items if row['kind']=='complete_switching_subfamily'))
        base.verify(a,record)
        partial=rotation_identity(a,np.asarray(record['s']),np.asarray(record['d']))
        transformed=dict(p=list(range(len(a))),s=[1]*len(a),
            d=(np.asarray(record['s'])*np.asarray(record['d'])).tolist(),cap=record['cap'])
        base.verify(partial,transformed)
        record.update(label=label,full_twisted_family_minimum=f['cap'],
            permutations_indispensable_for_full_minimum=record['cap']>f['cap'],
            equivalent_untwisted_partial_child=transformed,
            histogram=items[-1],events=items,
            source=str(source.relative_to(base.ROOT)),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
            included_source=str(dependency.relative_to(base.ROOT)),included_source_sha256=hashlib.sha256(dependency.read_bytes()).hexdigest(),
            evidence='all2^(n-1) switching signs and all2^n matching signs, with p=id; no claim beyond this exact subfamily')
        records.append(record)
        print(label,'switching',record['cap'],'full',f['cap'],flush=True)
    rng=np.random.default_rng(2026091902);identity_tests=[]
    for n in range(3,11):
        for k in range(15):
            upper=np.triu(rng.choice([-1,1],size=(n,n)),1);a=upper+upper.T
            s=rng.choice([-1,1],size=n);d=rng.choice([-1,1],size=n)
            rotation_identity(a,s,d)
            identity_tests.append(dict(n=n,A=a.tolist(),s=s.tolist(),d=d.tolist()))
    reversals=partial_reversal_caps(np.asarray(conference['child_matrix']))
    assert reversals['minimum_child_cap']>=15
    result=dict(records=records,random_signed_permutation_identity_tests=identity_tests,
        identity_test_count=len(identity_tests),identity_all_passed=True,
        identity='rotate pairs (x_i,y_i)=(yprime_i,-xprime_i) at s_i=-1; reverse only child edges internal to this set; matching becomes d*s',
        conference10_all_partial_reversals=reversals)
    path=base.ROOT/'computations/results/twisted_chiral_switching_audit_2026_09_19.json'
    path.write_text(json.dumps(result,indent=2)+'\n')
    print('conference partial reversal histogram',reversals['cap_histogram'])
    untwisted={row['label']:row['cap'] for row in json.loads((base.ROOT/'computations/results/twisted_chiral_untwisted_audit_2026_09_19.json').read_text())['records']}
    lines=['# Switching versus permutation in the finite twisted family','',
        'All matchings d are optimized in each column. Every fixed-child family comparison is exact; these are not global parent optima.','',
        '| Child | Untwisted | Switching only | Full family | Permutations indispensable for full minimum |',
        '|---|---|---|---|---|']
    for row in records:
        lines.append(f'| {row["label"]} | {untwisted[row["label"]]} | {row["cap"]} | {row["full_twisted_family_minimum"]} | {"yes" if row["permutations_indispensable_for_full_minimum"] else "no"} |')
    lines+=['','Switching-only means p=id and every2^(n-1) switching vector modulo global sign; every2^n matching d is tested.','',
        'At J={i:s_i=-1}, rotate (x_i,y_i)=(yprime_i,-xprime_i). The resulting child reverses only edges internal to J, the bridge becomes an untwisted copy of that child, and the matching becomes d*s. This exact signed-permutation matrix identity passed120 independently generated random examples.','',
        'All1024 partial reversals of the conference order10 child have Q histogram {15:82,17:762,19:180}; none is an optimal Q13 child.','',
        'Source/results: computations/twisted_chiral_switching_audit_2026_09_19.py and computations/results/twisted_chiral_switching_audit_2026_09_19.json.','']
    (base.ROOT/'artifacts/twisted_chiral_switching_table_2026_09_19.md').write_text('\n'.join(lines))


if __name__=='__main__':main()
