"""Exact H2 lift of stored small minimizers; no claims beyond tested matrices."""
import hashlib
import json
import math
from pathlib import Path
import subprocess
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
CPP=ROOT/'computations/decisive_independent_h2_lift_caps_2026_09_07.cpp'
BIN=CPP.with_suffix('')
subprocess.run(['g++','-O3','-std=c++17',str(CPP),'-o',str(BIN)],check=True)

def scan(a):
    bits=''.join('1' if a[i,j]>0 else '0' for i in range(len(a)) for j in range(i+1,len(a)))
    return json.loads(subprocess.check_output([str(BIN),str(len(a)),bits],text=True))

def fromcode(n,code):
    a=np.ones((n,n),dtype=int);np.fill_diagonal(a,0);k=0
    for i in range(1,n):
        for j in range(i+1,n):
            a[i,j]=a[j,i]=1 if ((code>>k)&1)==0 else -1;k+=1
    return a

inputs=[]
for n,code in [(3,0),(4,1),(5,13),(6,220),(7,828),(7,826),(8,53014),(9,898008),(9,6737136)]:
    inputs.append((f'gauged_n{n}_code{code}',fromcode(n,code)))
for n,code in [(5,13),(6,220)]:
    a=fromcode(n,code);d=np.eye(n,dtype=int)
    inputs.append((f'H2_lift_of_gauged_n{n}_code{code}',np.block([[a,a+d],[a+d,-a]])))
for name in ['exact_m10.json','nested_10_in_11_cap17.json','heuristic_m11.json','extension_nested_m11_to_12.json']:
    path=ROOT/'computations/results'/name;data=json.loads(path.read_text())
    key=next(k for k in ['matrix','parent_matrix','conference_matrix'] if k in data)
    a=np.array(data[key],dtype=int);inputs.append((str(path.relative_to(ROOT)),a))
    if 'nested_10' in name:inputs.append((name+' first10',a[:10,:10]))

records=[]
for name,a in inputs:
    n=len(a);base=scan(a);d=np.eye(n,dtype=int)
    b=np.block([[a,a+d],[a+d,-a]])
    result=scan(b)
    spin=np.ones(2*n,dtype=int)
    for j in range(1,2*n):spin[j]=1-2*((result['max_gray_spin_code']>>(j-1))&1)
    x,y=spin[:n],spin[n:]
    pp=(x+y)//2;rr=(x-y)//2
    uu=int(pp@a@pp)//2;vv=int(rr@a@rr)//2;ww=int(pp@a@rr)
    dd=int(x@y)
    assert 2*(uu-vv+ww)+dd==result['max']
    record=dict(source=name,parent=base,lift=result,
        lift_ratio=result['cap']/(2*n)**1.5,parent_ratio=base['cap']/n**1.5,
        delta_from_homogeneous=result['cap']-2*math.sqrt(2)*base['cap'],
        delta_per_parent_vertex=(result['cap']-2*math.sqrt(2)*base['cap'])/n,
        maximizing_support=dict(S_size=int(np.count_nonzero(pp)),T_size=int(np.count_nonzero(rr)),
            u=uu,v=vv,w=ww,diagonal_correction=dd,
            old_q_x=uu+vv+ww,old_q_y=uu+vv-ww),
        lift_matrix=b.tolist())
    records.append(record);print(json.dumps({k:v for k,v in record.items() if k!='lift_matrix'}),flush=True)
out=ROOT/'computations/decisive_independent_h2_lift_caps_2026_09_07.json'
out.write_text(json.dumps(dict(status='exact finite cap scans; ratios floating diagnostics',
    cpp_sha256=hashlib.sha256(CPP.read_bytes()).hexdigest(),records=records),indent=2)+'\n')
