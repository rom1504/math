import numpy as np,math,itertools
from scipy.optimize import differential_evolution,minimize
z=np.load('/home/math/quadra/tmp/rowcert/cand.npz');L=z['L'];lam=5.382104195764755;m=4;n=4
w=np.exp(-lam*(L-L.min()));q=w/w.sum();idx=np.arange(len(L),dtype=np.uint64);rows=np.stack([((idx>>(i*n))&15).astype(int) for i in range(m)])
features=[2,1,4,8]
ys=np.stack([(np.bitwise_count((rows[i]&features[i]).astype(np.uint8))%2).astype(int) for i in range(m)])
code=sum(ys[i]<<i for i in range(m));Q=np.bincount(code,weights=q,minlength=16);print('Q',Q,'min',Q.min(),'H',-np.dot(Q,np.log(Q)))
f=-np.log(Q).reshape((2,)*m,order='F') # because code i bits / axis i, reshape F
# expectation tensor product pi where pi_i=[1-a,a]
def fg(a):
 ps=[np.array([1-x,x]) for x in a];cur=f
 for i in reversed(range(m)):cur=np.tensordot(cur,ps[i],axes=(i,0))
 val=float(cur)-sum(-(x*math.log(x)+(1-x)*math.log(1-x)) if 0<x<1 else 0 for x in a)
 # easier numerical grad finite no
 return val
res=differential_evolution(fg,[(1e-12,1-1e-12)]*m,tol=1e-12,popsize=50,maxiter=3000,polish=True,seed=1)
print('best',res.fun,res.x,'point masses',-np.log(Q).min(),-np.log(Q).max())
# coordinate descent all starts
def cd(a):
 a=np.array(a,float)
 for z in range(10000):
  old=a.copy()
  for i in range(m):
   others=[np.array([1-a[j],a[j]]) for j in range(m)]
   eff=f
   for j in reversed(range(m)):
    if j==i:continue
    eff=np.tensordot(eff,others[j],axes=(j,0))
   a[i]=1/(1+math.exp(eff[1]-eff[0]))
  if max(abs(a-old))<1e-14:break
 return fg(a),a
vals=[]
for seed in range(100):
 rng=np.random.default_rng(seed);vals.append(cd(rng.random(m)))
print('cd distinct',sorted({round(x[0],12) for x in vals}),sorted(vals,key=lambda x:x[0])[:5])
# C matrix rectangle oscillation exact binary: max over rest differences
C=np.zeros((m,m))
for i in range(m):
 for j in range(i):
  rest=[k for k in range(m) if k not in [i,j]];best=0
  for bits in itertools.product([0,1],repeat=len(rest)):
   sl=[slice(None)]*m
   for k,b in zip(rest,bits):sl[k]=b
   M=f[tuple(sl)]
   # axes remain original order among i,j perhaps 2x2
   v=abs(M[0,0]-M[0,1]-M[1,0]+M[1,1]);best=max(best,v)
  C[i,j]=C[j,i]=best
print('C',C,'rho',np.linalg.eigvalsh(C)[-1],'condition rho<4',np.linalg.eigvalsh(C)[-1]<4)
