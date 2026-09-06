"""Restricted rank-one weave tests. Exact identities and lower witnesses only."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import numpy as np

ROOT=Path(__file__).resolve().parents[1]

def walsh(m):
    assert m>0 and m&(m-1)==0
    par=np.array([bin(i).count('1')%2 for i in range(m)],dtype=np.int64)
    return 1-2*par[np.bitwise_and(np.arange(m)[:,None],np.arange(m)[None,:])]

def instance(m,k,seed):
    rng=np.random.default_rng(seed)
    h=walsh(m)
    perms=np.stack([rng.permutation(m) for _ in range(m)])
    rows=np.stack([rng.choice(m,k,replace=False) for _ in range(m)])
    t=np.stack([h[rows[i]][:,perms[i]] for i in range(m)])
    s=rng.choice([-1,1],size=(m,m))
    s=np.triu(s)+np.triu(s,1).T
    b=np.einsum('iaj,jbi,ij->iajb',t,t,s).reshape(m*k,m*k)
    assert np.array_equal(b,b.T) and np.all(np.abs(b)==1)
    a=b.copy();np.fill_diagonal(a,0)
    initials=[]
    for _ in range(8):
        order=rng.permutation(m); pi=np.empty(m,dtype=int);u=np.ones(m,dtype=int)
        for j in range(0,m,2):
            v,w=order[j:j+2];pi[v]=w;pi[w]=v;u[w]=s[v,w]
        initials.append(np.concatenate([u[i]*t[i,:,pi[i]] for i in range(m)]))
    return a,b,t,s,perms,rows,initials

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--m',type=int,default=16)
    ap.add_argument('--k',type=int,default=14);ap.add_argument('--seed',type=int,default=9061001)
    ap.add_argument('--runs',type=int,default=200);ap.add_argument('--sweeps',type=int,default=80)
    ap.add_argument('--output',type=Path,required=True);args=ap.parse_args()
    m,k=args.m,args.k;n=m*k
    a,b,t,s,perms,rows,initials=instance(m,k,args.seed)
    binary=ROOT/'tmp/continued_convergence_weave_ascent_2026_09_06'
    source=ROOT/'computations/continued_convergence_weave_ascent_2026_09_06.cpp'
    subprocess.run(['g++','-O3','-std=c++17',str(source),'-o',str(binary)],check=True)
    data=f'{n} {args.runs} {args.sweeps} {args.seed}\n'
    data+='\n'.join(' '.join(map(str,r)) for r in a)+'\n'+str(len(initials))+'\n'
    data+='\n'.join(' '.join(map(str,r)) for r in initials)+'\n'
    result=json.loads(subprocess.run([str(binary)],input=data,text=True,capture_output=True,check=True).stdout)
    x=np.array(result['witness'],dtype=np.int64)
    exact=int(x@a@x)//2;assert abs(exact)==result['lower_Q']
    h=np.einsum('iaj,ia->ij',t,x.reshape(m,k))
    energy=int(np.sum(s*h*h.T));assert energy==int(x@b@x)
    norm=m*m*k;assert int(np.sum(h*h))==norm
    defects={}
    for sg in [-1,1]:
        defect=int(np.sum((h-sg*s*h.T)**2))
        assert defect==2*(norm-sg*energy)
        defects[str(sg)]=defect/(2*norm)
    result.update({'m':m,'k':k,'N':n,'p':k/m,'seed':args.seed,'runs':args.runs,
      'sweeps':args.sweeps,'normalized_lower_Q':result['lower_Q']/n**1.5,
      'full_rayleigh_abs':abs(energy)/norm,'defect_fractions':defects,
      'subhalf_required_defect':1-math.sqrt(k/m),'spectral_normalized_upper':.5*math.sqrt(m/k)+.5/math.sqrt(n),
      'row_selectors':rows.tolist(),'column_permutations':perms.tolist(),'S':s.tolist(),
      'status':'lower witness only; not an upper cap or an asymptotic claim'})
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({q:result[q] for q in ['m','k','N','normalized_lower_Q','full_rayleigh_abs','defect_fractions','subhalf_required_defect']}))

if __name__=='__main__':main()
