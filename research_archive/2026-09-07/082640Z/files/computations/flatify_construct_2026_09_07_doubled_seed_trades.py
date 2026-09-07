"""Actual doubled-seed flat sign matrices under orthogonal closed-quad trades."""
import argparse
import itertools
import json
from pathlib import Path
import numpy as np


def options(h):
    n = len(h)
    groups = {}
    for i,j in itertools.combinations(range(n),2):
        groups.setdefault(np.packbits(h[i] != h[j]).tobytes(), []).append((i,j))
    out = set()
    for pairs in groups.values():
        for ab,cd in itertools.combinations(pairs,2):
            q = tuple(sorted(ab+cd))
            if len(set(q)) == 4 and int(h[np.ix_(q,q)].sum()) % 8 == 0:
                out.add(q)
    return sorted(out)


def trade(h,q):
    q = list(q)
    r = np.ones((4,4), dtype=np.int64)-2*np.eye(4,dtype=np.int64)
    new = h.copy()
    new[q] = (r @ new[q])//2
    new[:,q] = (new[:,q] @ r)//2
    assert np.all(np.abs(new)==1) and np.array_equal(new,new.T)
    assert np.trace(new)==np.trace(h)
    return new


def score(h,x):
    e = np.einsum('bi,ij,bj->b',x,h,x,optimize=True)//2
    e -= np.trace(h)//2
    cap = int(np.max(np.abs(e)))
    count = int(np.count_nonzero(np.abs(e)==cap))
    return cap,count


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--m',type=int,default=8)
    parser.add_argument('--steps',type=int,default=100)
    parser.add_argument('--sample',type=int,default=12)
    parser.add_argument('--opt-diagonal',action='store_true')
    args=parser.parse_args()
    m=args.m;n=2*m
    src=Path(f'computations/results/exact_m{m}.json')
    a=np.array(json.loads(src.read_text())['matrix'],dtype=np.int64)
    s=a+np.eye(m,dtype=np.int64)
    h=np.kron(np.array([[1,1],[1,-1]],dtype=np.int64),s)
    x=1-2*((np.arange(1<<(n-1),dtype=np.int64)[:,None]>>np.arange(n))&1)
    selected_diagonal=np.ones(m,dtype=np.int64)
    if args.opt_diagonal:
        base=np.kron(np.array([[1,1],[1,-1]],dtype=np.int64),a)
        energies=np.einsum('bi,ij,bj->b',x,base,x,optimize=True)//2
        z=x[:,:m]*x[:,m:]
        ids=((z<0)*(1<<np.arange(m))).sum(axis=1)
        maxima=np.full(1<<m,-n*n,dtype=np.int64)
        minima=np.full(1<<m,n*n,dtype=np.int64)
        np.maximum.at(maxima,ids,energies)
        np.minimum.at(minima,ids,energies)
        all_d=1-2*((np.arange(1<<m)[:,None]>>np.arange(m))&1)
        linear=all_d@all_d.T
        caps=np.maximum(np.max(maxima[None,:]+linear,axis=1),
                        np.max(-minima[None,:]-linear,axis=1))
        selected_diagonal=all_d[int(np.argmin(caps))]
        s=a+np.diag(selected_diagonal)
        h=np.kron(np.array([[1,1],[1,-1]],dtype=np.int64),s)
        assert score(h,x)[0]==int(caps.min())
    initial=h.copy()
    rng=np.random.default_rng(20260907+m)
    best=score(h,x);records=[dict(step=0,score=best,matrix=h.tolist())]
    log=[]
    for step in range(args.steps+1):
        current=score(h,x)
        opts=options(h)
        row=dict(step=step,score=current,options=len(opts),
                 edits=int(np.count_nonzero(np.triu(h!=initial,1))))
        log.append(row)
        if step%10==0:print(json.dumps(row),flush=True)
        if step==args.steps or not opts:break
        if len(opts)>args.sample:
            opts=[opts[i] for i in rng.choice(len(opts),args.sample,replace=False)]
        candidates=[]
        for q in opts:
            new=trade(h,q)
            candidates.append((score(new,x),q,new))
        candidates.sort(key=lambda item:item[0])
        # Mostly cap-aware descent, occasional exploration across ties/barriers.
        chosen=candidates[0] if rng.random()<.85 else candidates[int(rng.integers(len(candidates)))]
        current,q,h=chosen
        if current<best:
            best=current
            records.append(dict(step=step+1,score=best,matrix=h.tolist()))
            print(json.dumps(dict(improved=True,step=step+1,score=best)),flush=True)
    seed_x=1-2*((np.arange(1<<(m-1))[:,None]>>np.arange(m))&1)
    seedcap=score(a,seed_x)[0]
    output=dict(m=m,n=n,seed_cap=seedcap,selected_diagonal=selected_diagonal.tolist(),
                aligned_child_target=2*np.sqrt((n-1)/(m-1))*seedcap,
                best=best,records=records,log=log)
    suffix='_optdiag' if args.opt_diagonal else ''
    Path(f'computations/results/flatify_construct_2026_09_07_doubled_seed_trades_m{m}{suffix}.json').write_text(
        json.dumps(output,indent=2)+'\n')


if __name__=='__main__':main()
