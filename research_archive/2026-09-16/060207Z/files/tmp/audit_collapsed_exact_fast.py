#!/usr/bin/env python3
"""Vectorized exact collapsed-shore recoupling defect by selected mask."""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

ROOT=Path('/home/math/quadra')
FILES=['conference_double_p5.json','conference_order10_gf9.json','conference_double_p13.json','conference_double_p17.json']

def spins(n):
 c=np.arange(1<<max(0,n-1),dtype=np.uint32)[:,None]
 b=((c>>np.arange(max(0,n-1),dtype=np.uint32))&1).astype(np.int16)
 return np.c_[np.ones(len(c),np.int16),1-2*b] if n else np.empty((1,0),np.int16)

def main():
 for fn in FILES:
  C=np.array(json.load(open(ROOT/'computations/results'/fn))['conference_matrix'],dtype=np.int16); n=len(C); XX=spins(n)
  groups={}
  for X in XX:
   ell=(C@X)*X; I=tuple(np.flatnonzero(ell>0)); groups.setdefault(I,[]).append((X,ell))
  total=0; maxdef=0; hist={}
  for I,rows in groups.items():
   I=np.array(I,dtype=int); J=np.array([a for a in range(n) if a not in set(I)],dtype=int); RR=spins(len(J)).astype(np.int32)
   for X,ell in rows:
    B=C*X[:,None]*X[None,:]; D=B[np.ix_(J,J)].astype(np.int32); h=B[np.ix_(J,I)].astype(np.int32)@np.ones(len(I),np.int32)
    p=int(np.ones(len(I),np.int32)@B[np.ix_(I,I)]@np.ones(len(I),np.int32)); r=int(np.ones(len(J),np.int32)@D@np.ones(len(J),np.int32)); V=p-r
    vals=np.einsum('bi,ij,bj->b',RR,D,RR,optimize=True)+2*np.abs(RR@h)
    aligned=int(vals.max()) # P positive in nontrivial row law; for degenerate this clips below
    defect=max(0,V-(p+aligned)); total+=defect; maxdef=max(maxdef,defect); hist[defect]=hist.get(defect,0)+1
  print(json.dumps({'file':fn,'n':n,'samples':len(XX),'unique_shores':len(groups),'mean_defect':total/len(XX),'mean_defect_over_n32':total/len(XX)/n**1.5,'max_defect':maxdef,'hist':hist},sort_keys=True),flush=True)
if __name__=='__main__':main()
