#!/usr/bin/env python3
"""Check explicit free-shore choices for the row-sign conference law."""

from __future__ import annotations

import json
from pathlib import Path
import numpy as np

ROOT=Path('/home/math/quadra')
FILES=['conference_double_p5.json','conference_order10_gf9.json','conference_double_p13.json','conference_double_p17.json']

def sg(v):
    return np.where(v>=0,1,-1).astype(np.int16)

def main():
  for fn in FILES:
    C=np.array(json.load(open(ROOT/'computations/results'/fn))['conference_matrix'],dtype=np.int16); n=len(C)
    codes=np.arange(1<<(n-1),dtype=np.uint32)[:,None]; bits=((codes>>np.arange(n-1,dtype=np.uint32))&1).astype(np.int16); XX=np.c_[np.ones(len(codes),np.int16),1-2*bits]
    sums={k:0.0 for k in ['response','P','minusR','loss_r1','loss_signh','loss_signDh','loss_signellJ','loss_signD1','best_explicit']}
    maxloss={k:0 for k in sums if k.startswith('loss')}
    for X in XX:
      ell=(C@X)*X; g=sg(ell); I=np.flatnonzero(g>0); J=np.flatnonzero(g<0); B=C*X[:,None]*X[None,:]
      D=B[np.ix_(J,J)].astype(np.int64); h=B[np.ix_(J,I)]@np.ones(len(I),np.int64)
      P=int(np.ones(len(I),np.int64)@B[np.ix_(I,I)]@np.ones(len(I),np.int64)); R=int(np.ones(len(J),np.int64)@D@np.ones(len(J),np.int64)); V=P-R
      vals=[]
      for name,r in [('r1',np.ones(len(J),np.int64)),('signh',sg(h)),('signDh',sg(D@h)),('signellJ',sg(ell[J])),('signD1',sg(D@np.ones(len(J),np.int64)))]:
        w=int(r@D@r+2*abs(h@r)); cert=P+w
        loss=max(0,V-cert); sums['loss_'+name]+=loss; maxloss['loss_'+name]=max(maxloss['loss_'+name],loss); vals.append(cert)
      sums['best_explicit']+=max(vals); sums['response']+=V; sums['P']+=P; sums['minusR']+=-R
    out={'file':fn,'n':n,'samples':len(XX),'means':{k:v/len(XX) for k,v in sums.items()},'maxloss':maxloss}
    print(json.dumps(out,sort_keys=True),flush=True)

if __name__=='__main__': main()
