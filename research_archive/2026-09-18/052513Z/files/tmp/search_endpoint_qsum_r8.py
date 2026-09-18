#!/usr/bin/env python3
"""MILP search for a violation of Q_S+Q_T <= R-|H|.

Gauge p=1 and prescribe the negative endpoint as the S|T cut.  Endpoint
optimality is exactly 0 <= C(R) <= c=C(T) for every cut R.  For fixed
candidate child cuts and signs, maximize the corresponding linear lower
bound |E_S|+|E_T|+|H|-4c.
"""

from itertools import combinations, product
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds


def edge_list(n): return [(i,j) for i in range(n) for j in range(i+1,n)]


def vec_sum(n, edge_pred):
    # Coefficients and constant for sum_{edges pred} a_e, a_e=2x_e-1.
    es=edge_list(n); z=np.zeros(len(es))
    ids=[k for k,e in enumerate(es) if edge_pred(*e)]
    z[ids]=2
    return z,-len(ids)


def add(u,cu,v,cv,alpha=1): return u+alpha*v,cu+alpha*cv


def solve(n,sn,kwS,kwT,signS,signT,signH):
    S=set(range(sn));T=set(range(sn,n)); W=set(range(kwS)); Z=set(range(sn,sn+kwT))
    cross=lambda i,j:(i in S)!=(j in S)
    internal=lambda i,j:(i in S)==(j in S)
    c,cc=vec_sum(n,cross)
    it,ci=vec_sum(n,internal)
    H,cH=2*it,2*ci
    def child_energy(X,W):
        h,ch=vec_sum(n,lambda i,j:i in X and j in X)
        h,ch=2*h,2*ch
        k,ck=vec_sum(n,lambda i,j:i in X and j in X and ((i in W)!=(j in W)))
        return h-4*k,ch-4*ck
    eS,ceS=child_energy(S,W);eT,ceT=child_energy(T,Z)
    # objective maximize signed child energies + signed H -4c
    obj=signS*eS+signT*eT+signH*H-4*c
    cobj=-obj
    const=signS*ceS+signT*ceT+signH*cH-4*cc
    rows=[];lo=[];hi=[]
    # all cuts modulo complement: last vertex fixed outside
    for mask in range(1<<(n-1)):
        R={i for i in range(n-1) if mask>>i&1}
        cut,kc=vec_sum(n,lambda i,j:(i in R)!=(j in R))
        # cut >=0 and cut-c <=0
        rows.append(cut);lo.append(-kc);hi.append(np.inf)
        rows.append(cut-c);lo.append(-np.inf);hi.append(-(kc-cc))
    # chosen signs genuinely lower-bound absolute values
    for vec,con,sgn in ((eS,ceS,signS),(eT,ceT,signT),(H,cH,signH)):
        rows.append(sgn*vec);lo.append(-sgn*con);hi.append(np.inf)
    ans=milp(cobj,integrality=np.ones(len(cobj)),bounds=Bounds(0,1),constraints=LinearConstraint(np.array(rows),np.array(lo),np.array(hi)),options={'time_limit':20})
    if not ans.success:return None
    value=round(-ans.fun+const)
    return value,tuple(round(x) for x in ans.x)


def main():
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('n',type=int);z=ap.parse_args();n=z.n
    best=(-10**9,None)
    for sn in range(1,n):
      # subset sizes suffice by within-shore permutation symmetry
      for kwS in range(sn+1):
       for kwT in range(n-sn+1):
        for signs in product((-1,1),repeat=3):
         ans=solve(n,sn,kwS,kwT,*signs)
         if ans and ans[0]>best[0]:
          best=(ans[0],(sn,kwS,kwT,signs,ans[1]));print('BEST',best,flush=True)
    print('FINAL',best)


if __name__=='__main__':main()
