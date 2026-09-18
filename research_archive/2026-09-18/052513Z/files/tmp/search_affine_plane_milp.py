import argparse
import itertools
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix

def solve(n, counts, w=1, timelimit=60):
    types=[]
    for a,c in enumerate(counts): types += [a]*c
    assert len(types)==n
    edges=[(i,j) for i in range(n) for j in range(i+1,n)]
    ne=len(edges); nz=4*n; mi=ne+nz
    # variables y edge binary, z_i,h binary, M continuous
    nv=mi+1; mid=mi
    q=np.ones((4,ne),dtype=int)
    for h in range(4):
      for e,(i,j) in enumerate(edges):
        char=types[i]^types[j]^(w if h else 0) # wrong: dot below; w is vector bits
        # q_h=(-1)^dot(types_i xor types_j xor wvec,hvec)
        q[h,e]=1 if (bin(char & h).count('1')%2)==0 else -1
    # reject non-rank2 edge words
    if len({tuple(row) for row in q})<4:return None
    rows=[]; cols=[]; dat=[]; lo=[]; hi=[]
    def add(coeff, low, high):
      rr=len(lo)
      for c,v in coeff.items():
        if v: rows.append(rr);cols.append(c);dat.append(v)
      lo.append(low);hi.append(high)
    const_edge=np.ones(ne,dtype=int)
    # score_h = 2 q_h*y - sum q_h equals M
    for h in range(4):
      coeff={e:2*q[h,e] for e in range(ne)}; coeff[mid]=-1
      rhs=int(q[h].sum())
      add(coeff,rhs,rhs)
    # all cut energies in [-M,M]
    for mask in range(1<<(n-1)):
      x=np.ones(n,dtype=int)
      for i in range(1,n):
        if (mask>>(i-1))&1:x[i]=-1
      ce=np.array([x[i]*x[j] for i,j in edges],dtype=int)
      # 2 ce*y-sum ce <= M
      coeff={e:2*ce[e] for e in range(ne)};coeff[mid]=-1
      add(coeff,-np.inf,int(ce.sum()))
      # 2 ce*y-sum ce >= -M => 2ce*y+M >= sum ce
      coeff={e:2*ce[e] for e in range(ne)};coeff[mid]=1
      add(coeff,int(ce.sum()),np.inf)
    # local r_i(h)>=0 and if z_i,h=1 then r_i(h)=0 via r<= (n-1)(1-z)
    inc=[[] for _ in range(n)]
    for e,(i,j) in enumerate(edges):inc[i].append(e);inc[j].append(e)
    for i in range(n):
      for h in range(4):
        es=inc[i]
        csum=sum(q[h,e] for e in es)
        coeff={e:2*q[h,e] for e in es}
        add(coeff,csum,np.inf) # 2qy >= sum q means r>=0
        zidx=ne+4*i+h
        coeff[zidx]=n-1
        add(coeff,-np.inf,csum+n-1) # r+(n-1)z <= n-1
      add({ne+4*i+h:1 for h in range(4)},1,np.inf)
    A=coo_matrix((dat,(rows,cols)),shape=(len(lo),nv)).tocsc()
    c=np.zeros(nv); c[mid]=1
    integrality=np.zeros(nv);integrality[:mi]=1
    lb=np.zeros(nv);ub=np.ones(nv);lb[mid]=0;ub[mid]=ne
    res=milp(c,integrality=integrality,bounds=Bounds(lb,ub),constraints=LinearConstraint(A,lo,hi),
             options={'time_limit':timelimit,'mip_rel_gap':0})
    print('n',n,'counts',counts,'w',w,'status',res.message,'fun',res.fun)
    if res.x is None:return None
    avec=np.array([1 if z>.5 else -1 for z in res.x[:ne]],dtype=int)
    M=int(round(res.x[mid]))
    zs=np.rint(res.x[ne:ne+nz]).astype(int).reshape(n,4)
    # exact verify
    scores=q@avec
    fields=np.zeros((4,n),int)
    for h in range(4):
      for e,(i,j) in enumerate(edges):fields[h,i]+=avec[e]*q[h,e];fields[h,j]+=avec[e]*q[h,e]
    ens=[]
    for mask in range(1<<(n-1)):
      x=np.ones(n,dtype=int)
      for i in range(1,n):
        if (mask>>(i-1))&1:x[i]=-1
      ens.append(sum(avec[e]*x[i]*x[j] for e,(i,j) in enumerate(edges)))
    print('M',M,'actual',max(map(abs,ens)),'scores',scores.tolist(),'cover',np.any(fields==0,axis=0).tolist())
    print('types',types)
    print('fields')
    print(fields)
    print('A')
    AA=np.zeros((n,n),int)
    for e,(i,j) in enumerate(edges):AA[i,j]=AA[j,i]=avec[e]
    print(AA.tolist())
    return AA,types,fields,M

if __name__=='__main__':
  ap=argparse.ArgumentParser();ap.add_argument('--n',type=int,default=13);ap.add_argument('--counts',default='4,3,3,3');ap.add_argument('--w',type=int,default=1);ap.add_argument('--time',type=int,default=60)
  a=ap.parse_args(); solve(a.n,tuple(map(int,a.counts.split(','))),a.w,a.time)
