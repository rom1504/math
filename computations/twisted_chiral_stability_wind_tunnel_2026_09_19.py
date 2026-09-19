#!/usr/bin/env python3
"""Small seeded wind-tunnel; exact per-matrix counts, uncertified sample means.

Each positive full-parent one-spin maximum is counted projectively, checking
ALL coordinates (including the fixed enumeration coordinate zero). Full even
order means odd local fields, so strict positivity is correct. This script is
not a cap optimizer or an expectation certificate.
"""
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import random
import subprocess
import time
import numpy as np
from twisted_chiral_summarize_2026_09_18 import collect, FAMILY_TAGS

ROOT=Path(__file__).resolve().parents[1]
RESULTS=ROOT/'computations/results'
SCRATCH=ROOT/'tmp/twisted_chiral_2026_09_18'
SOURCE=ROOT/'computations/twisted_chiral_stability_histogram_2026_09_19.cpp'
BINARY=SCRATCH/'stability_histogram_2026_09_19'
SEED=2026091904

def matrix_input(a):
    return str(len(a))+'\n'+'\n'.join(' '.join(map(str,row)) for row in a)+'\n'

def histogram(a):
    return json.loads(subprocess.check_output([str(BINARY)],input=matrix_input(a),text=True))

def direct_histogram(a):
    a=np.asarray(a,dtype=np.int16);n=len(a);states=1<<(n-1)
    positive=Counter();stable=Counter();zeros=0;emin=1000;emax=-1000
    for begin in range(0,states,8192):
        bits=np.arange(begin,min(begin+8192,states),dtype=np.uint64)
        x=np.ones((len(bits),n),dtype=np.int16)
        x[:,1:]=1-2*((bits[:,None]>>np.arange(n-1,dtype=np.uint64))&1).astype(np.int16)
        products=x*(x@a);energies=products.sum(axis=1)//2
        emin=min(emin,int(energies.min()));emax=max(emax,int(energies.max()))
        zeros+=int(np.count_nonzero(energies==0))
        positive.update(map(int,energies[energies>0]))
        stable.update(map(int,energies[(energies>0)&np.all(products>0,axis=1)]))
    return dict(order=n,projective_states=states,minimum=emin,maximum=emax,cap=max(-emin,emax),zero_count=zeros,
        positive_energy_histogram={str(k):v for k,v in sorted(positive.items())},
        positive_stable_histogram={str(k):v for k,v in sorted(stable.items())})

def thresholds(hist,q,beta,familymin):
    tests={'doubling_2sqrt2Q':lambda e:e*e>8*q*q,
        'diamond_4beta_over3':lambda e:3*e>4*beta,
        'known_family_minimum':lambda e:e>familymin}
    return {name:dict(raw_positive_count=sum(v for e,v in hist['positive_energy_histogram'].items() if test(int(e))),
        positive_stable_count=sum(v for e,v in hist['positive_stable_histogram'].items() if test(int(e)))) for name,test in tests.items()}

def main():
    start=time.monotonic();SCRATCH.mkdir(parents=True,exist_ok=True)
    command=['g++','-O3','-march=native','-std=c++17',str(SOURCE),'-o',str(BINARY)]
    if not BINARY.exists() or BINARY.stat().st_mtime<SOURCE.stat().st_mtime:subprocess.run(command,check=True)
    rng=random.Random(SEED);families=collect(FAMILY_TAGS,'complete_family');groups=[]
    for r,cls in [(6,0),(8,0),(8,1),(10,0),(10,1)]:
        row=families[r,cls];groups.append((f'optimal_order{r}_class{cls}',row,[row]))
    conference=json.loads((RESULTS/'twisted_chiral_2026_09_18_alternative_conference.json').read_text())['records'][0]
    joint=json.loads((RESULTS/'twisted_chiral_2026_09_19_joint_target38.json').read_text())
    known=[conference];seen={conference['matrix_sha256']}
    for row in joint['records']:
        if row.get('cap')==40 and row['matrix_sha256'] not in seen:
            assert row['child_matrix']==conference['child_matrix']
            seen.add(row['matrix_sha256']);known.append(row)
    groups.append(('conference_order10',conference,known))
    records=[];summaries=[];direct_checks=[]
    for label,template,known_rows in groups:
        a=np.asarray(template['child_matrix'],dtype=np.int16);r=len(a);q=template['child_cap'];familymin=template['cap']
        spins=1-2*((np.arange(1<<r,dtype=np.uint64)[:,None]>>np.arange(r,dtype=np.uint64))&1).astype(np.int16)
        beta=int(np.abs(spins@a).sum(axis=1).max())
        samples=[]
        for k in range(24):
            p=list(range(r));rng.shuffle(p);s=[1]+[rng.choice([-1,1]) for _ in range(r-1)];d=[rng.choice([-1,1]) for _ in range(r)]
            samples.append(('random',k,p,s,d))
        for k,row in enumerate(known_rows):samples.append(('known_minimum',k,row['p'],row['s'],row['d']))
        group_records=[]
        for kind,k,p,s,d in samples:
            b=a[np.ix_(p,p)]*np.asarray(s)[:,None]*np.asarray(s)[None,:]
            c=b+np.diag(d);parent=np.block([[a,c],[c,-a]]).astype(np.int16)
            h=histogram(parent.tolist());assert h['cap']>=familymin
            assert h['minimum']==-h['maximum']
            assert 2*sum(h['positive_energy_histogram'].values())+h['zero_count']==h['projective_states']
            assert int(str(h['maximum'])) in [int(e) for e in h['positive_stable_histogram']]
            if kind=='known_minimum':assert h['cap']==familymin
            if kind=='known_minimum' and k==0:
                assert direct_histogram(parent)==h
                direct_checks.append(label)
            row=dict(group=label,kind=kind,sample_index=k,child_order=r,child_cap=q,child_beta=beta,
                known_full_family_minimum=familymin,p=p,s=s,d=d,child_matrix=a.tolist(),parent_matrix=parent.tolist(),
                matrix_sha256=hashlib.sha256(parent.astype(np.int8).tobytes()).hexdigest(),
                counts=h,threshold_counts=thresholds(h,q,beta,familymin))
            records.append(row);group_records.append(row)
        random_rows=[row for row in group_records if row['kind']=='random'];aggregate={}
        for name in random_rows[0]['threshold_counts']:
            raw=[row['threshold_counts'][name]['raw_positive_count'] for row in random_rows]
            stable=[row['threshold_counts'][name]['positive_stable_count'] for row in random_rows]
            assert [v>0 for v in raw]==[v>0 for v in stable]
            aggregate[name]=dict(sample_count=len(raw),mean_raw_count=sum(raw)/len(raw),mean_stable_count=sum(stable)/len(stable),
                raw_count_range=[min(raw),max(raw)],stable_count_range=[min(stable),max(stable)],
                fraction_with_violation=sum(v>0 for v in raw)/len(raw),
                ratio_total_stable_to_total_raw=sum(stable)/sum(raw) if sum(raw) else None)
        summary=dict(group=label,child_cap=q,child_beta=beta,familyminimum=familymin,threshold_aggregates=aggregate)
        summaries.append(summary);print(json.dumps(summary),flush=True)
        assert time.monotonic()-start<300,'authorized five-minute computation ceiling exceeded'
    result=dict(schema='twisted-chiral-stability-wind-tunnel-v1',generated_utc=datetime.now(timezone.utc).isoformat(),
        random_seed=SEED,random_sampling='24 independent uniform permutations, root-gauge uniform switching signs, uniform matching signs per group; known minima excluded from averages',
        quantifier='Each per-matrix histogram is exact. Sample averages are heuristic Monte Carlo diagnostics, NOT population expectations or upper-bound certificates.',
        thresholds='Strict positive energy above threshold: e^2>8Q(child)^2; 3e>4beta(child); e>known exact fixed-child family minimum.',
        stability='All parent coordinates, including projectively fixed coordinate0, require x_i(Dx)_i>0. Fields are odd and never zero.',
        direct_numpy_replay_groups=direct_checks,source_sha256=hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
        build_command=command,elapsed_seconds=time.monotonic()-start,summaries=summaries,records=records)
    output=RESULTS/'twisted_chiral_stability_wind_tunnel_2026_09_19.json'
    output.write_text(json.dumps(result,indent=2)+'\n');print(output,flush=True)

if __name__=='__main__':main()
