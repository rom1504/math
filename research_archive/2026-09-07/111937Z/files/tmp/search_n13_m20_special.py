import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import coo_matrix

n=13;m=6;X=list(range(6));z=6;Y=list(range(7,13));M=20
edges=[(i,j) for i in range(n) for j in range(i+1,n)]; ei={e:k for k,e in enumerate(edges)}
nv=len(edges)+12
rows=[];cols=[];dat=[];lo=[];hi=[]
def add(co,low,high):
 r=len(lo)
 for k,v in co.items():
  if v:rows.append(r);cols.append(k);dat.append(v)
 lo.append(low);hi.append(high)
def lin_signed(es, signs=None):
 # sum a_e * sign = 2 sum sign*y-sum sign, returns coeff,const subtracted form
 if signs is None:signs=[1]*len(es)
 return {ei[tuple(sorted(e))]:2*s for e,s in zip(es,signs)}, -sum(signs)

# Fix balanced singleton stars canonically.
sx=[1,1,1,-1,-1,-1]; sy=[1,1,1,-1,-1,-1]
for ii,s in zip(X,sx): add({ei[(ii,z)]:1},1 if s==1 else 0,1 if s==1 else 0)
for jj,s in zip(Y,sy): add({ei[(z,jj)]:1},1 if s==1 else 0,1 if s==1 else 0)

# Plane Fourier endpoint constraints: cross total M and internal totals cancel.
cross=[(i,j) for i in X for j in Y]
co,c0=lin_signed(cross); add(co,M-c0,M-c0)
ints=list(__import__('itertools').combinations(X,2))+list(__import__('itertools').combinations(Y,2))
co,c0=lin_signed(ints);add(co,-c0,-c0)

# b_i=|a_i|+1 and c_j=|d_j|+1. Binary u selects nonnegative internal row.
for side,other,off in [(X,Y,0),(Y,X,6)]:
 for kk,i in enumerate(side):
  ies=[tuple(sorted((i,j))) for j in side if j!=i]
  ces=[tuple(sorted((i,j))) for j in other]
  # Lplus=b-a-1 and Lminus=b+a-1, represented coeff + constant.
  cb,b0=lin_signed(ces); ca,a0=lin_signed(ies)
  cp=cb.copy()
  for e,v in ca.items():cp[e]=cp.get(e,0)-v
  p0=b0-a0-1
  cm=cb.copy()
  for e,v in ca.items():cm[e]=cm.get(e,0)+v
  m0=b0+a0-1
  u=len(edges)+off+kk; BIG=24
  # u=1 => Lplus=0 and a>=0; u=0 => Lminus=0 and a<=0
  # |Lplus| <= BIG(1-u)
  cp1=cp.copy();cp1[u]=BIG;add(cp1,-np.inf,BIG-p0)
  cp1=cp.copy();cp1[u]=-BIG;add(cp1,-BIG-p0,np.inf)
  # |Lminus| <= BIG*u
  cm1=cm.copy();cm1[u]=-BIG;add(cm1,-np.inf,-m0)
  cm1=cm.copy();cm1[u]=BIG;add(cm1,-m0,np.inf)
  # a>=-BIG(1-u): a-BIG*u >= -BIG
  aa=ca.copy();aa[u]=-BIG;add(aa,-BIG-a0,np.inf)
  # a<=BIG*u: a-BIG*u <=0
  aa=ca.copy();aa[u]=-BIG;add(aa,-np.inf,-a0)

# Exact absolute cap for every spin (mod global sign).
for mask in range(1<<(n-1)):
 x=[1]+[-1 if mask>>(i-1)&1 else 1 for i in range(1,n)]
 signs=[x[i]*x[j] for i,j in edges]
 co,c0=lin_signed(edges,signs)
 add(co,-M-c0,M-c0)

A=coo_matrix((dat,(rows,cols)),shape=(len(lo),nv)).tocsc()
res=milp(np.zeros(nv),integrality=np.ones(nv),bounds=Bounds(np.zeros(nv),np.ones(nv)),
 constraints=LinearConstraint(A,lo,hi), options={'time_limit':600,'mip_rel_gap':0,'presolve':True})
print(res.message)
if res.x is not None:
 aa=np.zeros((n,n),int)
 for k,(i,j) in enumerate(edges):aa[i,j]=aa[j,i]=1 if res.x[k]>.5 else -1
 print(aa.tolist())
