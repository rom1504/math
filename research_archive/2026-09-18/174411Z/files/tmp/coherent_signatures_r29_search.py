#!/usr/bin/env python3
"""Exact searches for local common-coset obstructions in A6/A8/A9."""

from __future__ import annotations

import itertools
import numpy as np

A6 = np.array([
 [0,1,1,1,1,1], [1,0,-1,-1,1,1], [1,-1,0,1,-1,1],
 [1,-1,1,0,1,-1], [1,1,-1,1,0,-1], [1,1,1,-1,-1,0]], dtype=np.int64)
A8 = np.array([
 [0,1,1,1,1,1,1,1], [1,0,1,-1,1,1,-1,-1],
 [1,1,0,1,-1,1,-1,-1], [1,-1,1,0,-1,-1,-1,1],
 [1,1,-1,-1,0,-1,1,-1], [1,1,1,-1,-1,0,1,1],
 [1,-1,-1,-1,1,1,0,1], [1,-1,-1,1,-1,1,1,0]], dtype=np.int64)
A9 = np.array([
 [0,1,-1,1,1,1,1,-1,-1], [1,0,1,-1,-1,-1,1,-1,-1],
 [-1,1,0,1,1,1,1,1,-1], [1,-1,1,0,1,1,1,-1,1],
 [1,-1,1,1,0,1,-1,1,-1], [1,-1,1,1,1,0,-1,-1,-1],
 [1,1,1,1,-1,-1,0,1,1], [-1,-1,1,-1,1,-1,1,0,-1],
 [-1,-1,-1,1,-1,-1,1,-1,0]], dtype=np.int64)


def set_partitions_with_sizes(n: int, sizes: tuple[int, ...]):
    """Generate each unlabeled partition with the prescribed block sizes."""
    sizes=tuple(sorted(sizes))
    def rec(remaining, todo, blocks):
      if not todo:
        yield tuple(tuple(B) for B in blocks); return
      size=todo[0]
      equal_next=len(todo)>1 and todo[1]==size
      anchor=remaining[0] if equal_next else None
      pool=remaining[1:] if anchor is not None else remaining
      choose=size-1 if anchor is not None else size
      for rest in itertools.combinations(pool,choose):
        block=((anchor,)+rest) if anchor is not None else rest
        block_set=set(block)
        nxt=tuple(i for i in remaining if i not in block_set)
        yield from rec(nxt,todo[1:],blocks+[tuple(sorted(block))])
    yield from rec(tuple(range(n)),sizes,[])

def x_of(mask: int, n: int) -> np.ndarray:
    # coordinate 0 is fixed +1; bit i-1 means coordinate i is -1.
    return np.array([1] + [-1 if (mask >> (i-1)) & 1 else 1 for i in range(1,n)], dtype=np.int64)

def canonical_mask(x: np.ndarray) -> int:
    if x[0] < 0:
        x = -x
    return sum((int(x[i]) < 0) << (i-1) for i in range(1,len(x)))

def product_mask(a: int, b: int, n: int) -> int:
    return canonical_mask(x_of(a,n)*x_of(b,n))

def data(A: np.ndarray, m: int, cap: int):
    n=len(A); masks=list(range(1<<(n-1))); xs={a:x_of(a,n) for a in masks}
    rows={a:int((A@xs[a])@(A@xs[a])) for a in masks}
    good={a for a in masks if rows[a] <= cap}
    cand={}
    for S in itertools.combinations(range(n),m):
        idx=np.array(S); B=A[np.ix_(idx,idx)]
        energies={a:abs(int(xs[a][idx]@B@xs[a][idx])) for a in masks}
        q=max(energies.values())
        cand[S]={a for a in good if energies[a]==q}
    return rows,good,cand

def triangle_search(A: np.ndarray, name: str, m: int, cap: int):
    n=len(A); rows,good,cand=data(A,m,cap)
    total=eligible=bad=0; first=None
    for core in itertools.combinations(range(n),m-1):
      outside=[i for i in range(n) if i not in core]
      for adds in itertools.combinations(outside,3):
        sels=[tuple(sorted(core+(a,))) for a in adds]; total+=1
        cs=[cand[S] for S in sels]
        if any(not C for C in cs): continue
        eligible+=1; exists=False
        pair={product_mask(a,b,n) for a in cs[0] for b in cs[1]}
        for p in pair:
          if any(product_mask(p,c,n) in good for c in cs[2]):
            exists=True; break
        if not exists:
          bad+=1
          if first is None: first=(sels,[len(C) for C in cs])
    print(name,m,cap,{"triangles":total,"eligible":eligible,"ternary_obstructions":bad,"first":first,
                       "good_spins":len(good),"candidate_range":(min(map(len,cand.values())),max(map(len,cand.values())))})
    if first is not None:
      sels,_=first; cs=[sorted(cand[S]) for S in sels]
      def word(a): return ''.join('+' if v==1 else '-' for v in x_of(a,n))
      detail=[]
      for S,C in zip(sels,cs):
        idx=np.array(S); B=A[np.ix_(idx,idx)]
        detail.append((S,[(word(a),rows[a],abs(int(x_of(a,n)[idx]@B@x_of(a,n)[idx]))) for a in C]))
      triples=[]
      for a in cs[0]:
       for b in cs[1]:
        for c in cs[2]:
         d=product_mask(product_mask(a,b,n),c,n)
         triples.append((word(a),word(b),word(c),word(d),rows[d]))
      print(' FIRST_DETAIL',detail)
      print(' TERNARY_PRODUCTS',triples)

def main():
    for A,name,settings in [
      (A6,"A6",[(3,30),(4,30),(5,30)]),
      (A8,"A8",[(3,40),(4,40),(4,64),(5,40),(5,64),(6,64),(7,64)]),
      (A9,"A9",[(4,40),(5,40),(5,80),(6,40),(6,80),(7,80),(8,80)])]:
      for m,cap in settings: triangle_search(A,name,m,cap)


def balanced_coset_helly(A: np.ndarray, m: int, cap: int, sizes=(1,2,2,2,2)):
    """Search Johnson triangles pair-covered but not jointly covered by row-good cosets."""
    n=len(A); rows,good,cand=data(A,m,cap)
    selectors=list(cand); si={S:i for i,S in enumerate(selectors)}
    # Here cand[S] already is exact-ground intersect row-good; a row-good coset
    # covers S iff it intersects cand[S].
    covers=[]; seen=set()
    for part in set_partitions_with_sizes(n,tuple(sizes)):
      reps=[B[0] for B in part]
      free=[i for i in range(n) if i not in reps]
      for bits in range(1<<len(free)):
        base=np.ones(n,dtype=np.int64)
        for j,i in enumerate(free):
          if (bits>>j)&1: base[i]=-1
        W=set()
        for zz in itertools.product((-1,1),repeat=len(part)):
          x=base.copy()
          for z,B in zip(zz,part): x[list(B)]*=z
          W.add(canonical_mask(x))
        key=tuple(sorted(W))
        if key in seen: continue
        seen.add(key)
        if not W <= good: continue
        cov=frozenset(i for i,S in enumerate(selectors) if W & cand[S])
        covers.append((key,cov,part))
    pair=set(); triple=set()
    for _,cov,_ in covers:
      for ab in itertools.combinations(sorted(cov),2): pair.add(ab)
      for abc in itertools.combinations(sorted(cov),3): triple.add(abc)
    total_star=eligible=walls=0; first=None
    for core in itertools.combinations(range(n),m-1):
      outside=[i for i in range(n) if i not in core]
      for adds in itertools.combinations(outside,3):
        total_star+=1
        ids=tuple(sorted(si[tuple(sorted(core+(a,)))] for a in adds))
        if all(tuple(sorted(e)) in pair for e in itertools.combinations(ids,2)):
          eligible+=1
          if ids not in triple:
            walls+=1
            if first is None: first=tuple(selectors[i] for i in ids)
    result={'n':n,'m':m,'cap':cap,'sizes':sizes,'unique_cosets':len(seen),
      'row_good_cosets':len(covers),'pairwise_triangles':eligible,'walls':walls,'first':first}
    result['star_triangles']=total_star
    print('BALANCED_HELLY',result)
    if first is not None:
      ids=tuple(si[S] for S in first)
      def word(a): return ''.join('+' if v==1 else '-' for v in x_of(a,n))
      for uv in itertools.combinations(ids,2):
        key,cov,part=next(record for record in covers if set(uv)<=record[1])
        witnesses=[]
        for i in uv:
          hits=sorted(set(key)&cand[selectors[i]])
          witnesses.append((selectors[i],[(word(a),rows[a]) for a in hits]))
        print(' EDGE_COSET',tuple(selectors[i] for i in uv),{'partition':part,
              'base':word(key[0]),'max_row':max(rows[a] for a in key),'row_hist':sorted({rows[a] for a in key}),
              'witnesses':witnesses})
    return result


def verify_target_wall():
    """Verify the exact A9 pairwise-but-not-global balanced-coset wall."""
    result=balanced_coset_helly(A9,6,80)
    target=((0,1,2,3,4,5),(0,1,3,4,5,6),(0,1,3,4,5,7))
    assert result['unique_cosets']==15120
    assert result['row_good_cosets']==452
    assert result['star_triangles']==504
    assert result['pairwise_triangles']==122
    assert result['walls']==6 and result['first']==target

    # Among the cap-80 exact-child candidate triples, the least possible full
    # row cap after a max-size-two, five-block refinement is exactly 88.
    rows,good,cand=data(A9,6,80)
    assert sorted(set(rows.values()))==[16,24,32,40,48,56,64,72,80,88,96,104,112,120,128]
    choices=[sorted(cand[S]) for S in target]
    assert list(map(len,choices))==[2,2,2]
    parts=list(set_partitions_with_sizes(9,(1,2,2,2,2)))
    optimum=10**9
    specific_optimum=10**9
    for aa in itertools.product(*choices):
      xs=[x_of(a,9) for a in aa]; base=xs[0]
      classes={}
      for i in range(9):
        classes.setdefault(tuple(int(x[i]*base[i]) for x in xs[1:]),[]).append(i)
      for part in parts:
        if not all(any(set(block)<=set(C) for C in classes.values()) for block in part):
          continue
        maximum=0
        for zz in itertools.product((-1,1),repeat=5):
          x=base.copy()
          for z,block in zip(zz,part): x[list(block)]*=z
          maximum=max(maximum,rows[canonical_mask(x)])
        optimum=min(optimum,maximum)
        if aa==(1,99,99): specific_optimum=min(specific_optimum,maximum)
    assert optimum==88
    assert specific_optimum==88
    special=[x_of(a,9) for a in (1,99,99)]
    special_classes={}
    for i in range(9):
      special_classes.setdefault(tuple(int(x[i]*special[0][i]) for x in special[1:]),[]).append(i)
    assert sorted(map(len,special_classes.values()))==[3,6]
    print('MIN_COMMON_BALANCED_CAP',optimum)

    # Exact signature-profile identity on the Johnson membership map: adjacent
    # selectors move by Hamming distance two, but all n coordinate signatures
    # are distinct.
    n,m=9,6
    sels=list(itertools.combinations(range(n),m))
    family=[np.array([1 if i in S else -1 for i in range(n)],dtype=np.int64) for S in sels]
    base=family[0]
    signatures={tuple(int(x[i]*base[i]) for x in family[1:]) for i in range(n)}
    assert len(signatures)==n
    for S,T in itertools.combinations(range(len(sels)),2):
      if len(set(sels[S])^set(sels[T]))==2:
        assert int(np.sum(family[S]!=family[T]))==2
    print('PASS coherent_signatures_r29_search')

if __name__ == '__main__': verify_target_wall()
