import numpy as np,itertools,time
z=np.load('/home/math/quadra/tmp/rowcert/cand.npz');T=z['T'];p=z['p'];m=4;q=16

def marg(S):
 cur=T; axes=list(range(m))
 for ax in reversed(range(m)):
  if ax in S:continue
  pos=axes.index(ax);cur=np.tensordot(cur,p[ax],axes=(pos,0));axes.pop(pos)
 # axes remain S in sorted order
 return cur
pairs=np.array(list(itertools.combinations(range(q),2)),int)
def krect(A):
 k=A.ndim
 if k==1:return float(A.max()-A.min())
 # chunk first pairs
 best=0
 for lo in range(0,len(pairs),8):
  pp=pairs[lo:lo+8]
  D=A[pp[:,0]]-A[pp[:,1]]
  # D shape chunk,q,...
  for axis in range(1,k-1):
   # difference along first q axis after pair dims, always axis=1
   D=D[:,pairs[:,0],...]-D[:,pairs[:,1],...]
   # this results prior pair dimensions accumulating after first? test shapes
  # now shape chunks,P...,q; range last
  val=np.max(D,axis=-1)-np.min(D,axis=-1)
  best=max(best,float(np.max(val)))
 return best
for k in range(2,5):
 for S in itertools.combinations(range(m),k):
  A=marg(S);st=time.time();K=krect(A);print(S,A.shape,K,'time',time.time()-st)
